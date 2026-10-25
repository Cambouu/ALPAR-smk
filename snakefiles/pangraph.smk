from pathlib import Path

# PanGraph experiment: graph features from the assemblies themselves, compared
# with snippy's SNPs. PanGraph 1.4.0 builds one pangenome graph of every strain
# plus the reference (no reads, no reference backbone, no labels), and
# pangraph_features.py turns the graph into ALPAR's long feature table: block
# presence plus SNPs, insertions and deletions, named by reference coordinates
# where the reference has one copy of the block.
#
# Each set in config pangraph_feature_sets (default SNP + GPA and GPA + PanGraph)
# is merged, pivoted once for all antibiotics and trained with the baseline
# combined_ml; pangraph_report compares the sets. The SNP and GPA tables are the
# baseline's. Not part of automatix; target experiment_pangraph. Every rule here
# writes under pangraph/ or experiments/pangraph/.

PANGRAPH_OUT_DIR = OUT_DIR / "pangraph"
PANGRAPH_LOGS_DIR = PANGRAPH_OUT_DIR / "logs"
EXPERIMENT_PANGRAPH_DIR = OUT_DIR / "experiments" / "pangraph"

PREPARE_PANGRAPH_INPUT = SCRIPTS_DIR / "prepare_pangraph_input.py"
PANGRAPH_FEATURES = SCRIPTS_DIR / "pangraph_features.py"

# alpar-smk-pangraph:1.0.0: the pangraph binary from neherlab/pangraph:1.4.0 plus
# the Snakemake worker (pangraph/containers/pangraph.Dockerfile in the thesis
# repository). PanGraph has no Conda package, so the build rule has no conda env.
PANGRAPH_IMAGE = "pangraph@sha256:a65d99d6fdf0166f9e2fcfba5327982b7d5c07141953484508a6b5c585ed04b0"

PANGRAPH_FEATURE_SETS = config.get("pangraph_feature_sets", {
    "snp_gpa": ["snp", "gpa"],
    "gpa_pangraph": ["gpa", "pangraph"],
})

wildcard_constraints:
    pangraph_set="(?:" + "|".join(re.escape(name) for name in PANGRAPH_FEATURE_SETS) + ")"

# Same model settings as the baseline target ml.
PANGRAPH_ML_SETTINGS = dict(
    random_seed=[42],
    test_size=[0.2],
    resampling_strategy=["cv"],
    model_type=["xgb"],
    feature_importance_analysis_strategy=["gini"],
)


# -----------------------
# Graph and Features
# -----------------------

rule pangraph_input:
    input:
        store = rules.rename_files.output.store,
        # Every strain of the run: the graph is built without labels.
        samples = lambda wc: expand(
            Path(rules.rename_files.output.store) / "{sample}",
            sample = get_sample_names(wc)
        ),
        # MG1655 as GenBank, the reference snippy uses, so positions share its
        # coordinates. FASTA_FILE is the Prokka protein database, not a genome.
        reference = GBFF_FILE,
    output:
        fasta = PANGRAPH_OUT_DIR / "input.fa.gz",
        samples = PANGRAPH_OUT_DIR / "input_samples.txt",
    log: PANGRAPH_LOGS_DIR / "pangraph_input.log"
    benchmark: BENCHMARKS_DIR / "pangraph_input.tsv"
    params:
        ref_name = lookup(dpath="pangraph_ref_name", within=config, default="MG1655"),
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        # The script sorts the assemblies by name: PanGraph's graph depends on input order.
        r"""
        printf '%s\n' {input.samples:q} > {output.samples:q}
        python {PREPARE_PANGRAPH_INPUT:q} --reference {input.reference:q} --ref-name {params.ref_name:q} \
            --out {output.fasta:q} @{output.samples:q} > {log:q} 2>&1
        """

rule pangraph_build:
    input: rules.pangraph_input.output.fasta
    output: PANGRAPH_OUT_DIR / "graph.json.gz"
    log: PANGRAPH_LOGS_DIR / "pangraph_build.log"
    benchmark: BENCHMARKS_DIR / "pangraph_build.tsv"
    container: CONTAINERS.format(PANGRAPH_IMAGE)
    # Locally, 10 threads kept only about 2.8 cores busy and 4 were 7% slower.
    threads: 8
    resources:
        # Extrapolated from 5-21 local genomes (RAM ~ N^0.81): about 126 GB for
        # 1,868 genomes, and wall time about 3 to 6 days.
        mem_mb = 256000,
        htcondor_request_mem_mb = 244141,
    shell:
        # --verify checks that every input sequence is rebuilt exactly from the graph.
        r"""
        pangraph build {input:q} -o {output:q} -j {threads} --verify --no-progress-bar > {log:q} 2>&1
        """

rule pangraph_features:
    input: rules.pangraph_build.output[0]
    output:
        features = PANGRAPH_OUT_DIR / "pangraph_features.tsv",
        carriers = PANGRAPH_OUT_DIR / "feature_carriers.tsv",
    log: PANGRAPH_LOGS_DIR / "pangraph_features.log"
    benchmark: BENCHMARKS_DIR / "pangraph_features.tsv"
    params:
        ref_name = lookup(dpath="pangraph_ref_name", within=config, default="MG1655"),
        # Label-free: carriers are counted over all strains. 1 keeps every
        # feature, like the snippy SNP table.
        min_carriers = lookup(dpath="pangraph_min_carriers", within=config, default=1),
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    resources:
        # The whole graph JSON is held in memory: about 44 GB extrapolated for 1,868 genomes.
        mem_mb = 96000,
        htcondor_request_mem_mb = 91553,
    shell:
        r"""
        python {PANGRAPH_FEATURES:q} {input:q} {output.features:q} --ref {params.ref_name:q} \
            --min-carriers {params.min_carriers} --carriers {output.carriers:q} > {log:q} 2>&1
        """

rule pangraph:
    localrule: True
    input: rules.pangraph_features.output
    output: touch(OUT_DIR / "flags" / "pangraph.done")


# -----------------------
# Feature Table and ML
# -----------------------

PANGRAPH_FAMILY_TABLES = {
    "snp": rules.binary_mutation_table.output[0],
    "gpa": rules.binary_gpa.output[0],
    "pangraph": rules.pangraph_features.output.features,
}

# All families are built from every strain, so one table serves every antibiotic;
# combined_ml takes each antibiotic's strains from its split.
rule merge_pangraph_set:
    localrule: True
    input: lambda wc: [PANGRAPH_FAMILY_TABLES[family] for family in PANGRAPH_FEATURE_SETS[wc.pangraph_set]]
    output: temp(EXPERIMENT_PANGRAPH_DIR / "{pangraph_set}" / "feature_table" / "merged_table.tsv"),
    threads: 1,
    shell:
        r"""
        cat {input:q} > {output:q}
        """

use rule pivot_merged_features_miller as pivot_pangraph_set with:
    group: "pivot_pangraph_set_batch"
    input: rules.merge_pangraph_set.output,
    output: EXPERIMENT_PANGRAPH_DIR / "{pangraph_set}" / "feature_table" / "merged_table_pivot.tsv",
    log: EXPERIMENT_PANGRAPH_DIR / "{pangraph_set}" / "feature_table" / "logs" / "pivot_pangraph_set.log",
    benchmark: BENCHMARKS_DIR / "pivot_pangraph_set_{pangraph_set}.tsv",


PANGRAPH_ML_STEM = str(EXPERIMENT_PANGRAPH_DIR / "{pangraph_set}" / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}")

use rule combined_ml as combined_ml_pangraph with:
    group: "combined_ml_pangraph_batch"
    input:
        binary_mutation_table = rules.pivot_pangraph_set.output[0],
        phenotype_table = rules.phenotype_dataframe_creator.output[0],
        train = lambda wildcards: expand(rules.split_train_test.output, split_category=["train"], antibiotic=wildcards.antibiotic)[0],
        test = lambda wildcards: expand(rules.split_train_test.output, split_category=["test"], antibiotic=wildcards.antibiotic)[0],
    output:
        best_params = PANGRAPH_ML_STEM + "_best_params.txt",
        model_file  = PANGRAPH_ML_STEM + "_model.sav",
        result      = PANGRAPH_ML_STEM + "_Result.txt",
        fia         = PANGRAPH_ML_STEM + ".txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_pangraph_{pangraph_set}_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.tsv"
    log: PANGRAPH_ML_STEM + ".log"
    params:
        feature_importance_analysis = True,
        save_model = True,
        param_grid_size = lookup(dpath="ml_param_grid_size", within=config, default="small"),
        # The baseline's table size always selects the low-memory search (on the
        # test split); the SNP + GPA table is close to the size threshold, so
        # every set is forced into the same search.
        param_grid_low_memory_mode = True,
        parallel_fits = lookup(dpath="ml_parallel_fits", within=config, default=1),


# -----------------------
# Report
# -----------------------

def pangraph_ml_outputs(name):
    return expand(
        getattr(rules.combined_ml_pangraph.output, name),
        pangraph_set=list(PANGRAPH_FEATURE_SETS),
        antibiotic=ANTIBIOTICS,
        **PANGRAPH_ML_SETTINGS,
    )

rule pangraph_report:
    input:
        results = pangraph_ml_outputs("result"),
        best_params = pangraph_ml_outputs("best_params"),
        fia = pangraph_ml_outputs("fia"),
        pivots = expand(rules.pivot_pangraph_set.output[0], pangraph_set=list(PANGRAPH_FEATURE_SETS)),
        carriers = rules.pangraph_features.output.carriers,
        samples = rules.pangraph_input.output.samples,
    output:
        markdown = EXPERIMENT_PANGRAPH_DIR / "report.md",
        table = EXPERIMENT_PANGRAPH_DIR / "report.tsv",
    log: PANGRAPH_LOGS_DIR / "pangraph_report.log"
    params:
        feature_sets = PANGRAPH_FEATURE_SETS,
        reference_set = lookup(dpath="pangraph_reference_set", within=config, default=next(iter(PANGRAPH_FEATURE_SETS))),
        min_carriers = lookup(dpath="pangraph_min_carriers", within=config, default=1),
        features_log = rules.pangraph_features.log[0],
        # Earlier runs shown for reference when their Result files exist; not rerun.
        reference_ml_dirs = {
            "SNP + GPA + panproteome (baseline ml/)": str(OUT_DIR / "ml"),
            "SNP + GPA (experiment/feature-sets)": str(OUT_DIR / "experiments" / "snp_gpa" / "ml"),
        },
        benchmarks_dir = str(BENCHMARKS_DIR),
        # Quinolone resistance genes, MG1655 (NC_000913.3) coordinates, both on the minus strand.
        loci = {"gyrA": [2336793, 2339420], "parC": [3163715, 3165973]},
        top_features = 20,
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    resources:
        mem_mb = 8000,
    script:
        SCRIPTS_DIR / "pangraph_report.py"

rule experiment_pangraph:
    localrule: True
    input: rules.pangraph_report.output
    output: touch(OUT_DIR / "flags" / "experiment_pangraph.done")

import math
from pathlib import Path

# PanGraph experiment: graph features from the assemblies themselves, compared
# with snippy's SNPs. PanGraph 1.4.0 builds one pangenome graph of every strain
# plus the reference (no reads, no reference backbone, no labels), from
# subgraphs of strain groups joined by `pangraph merge`, and
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
# repository). PanGraph has no Conda package, so its rules have no conda env.
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

# One build of all 1,868 genomes went over 244 GB within 32 minutes (job 68286,
# 2026-10-26). The strains are therefore split into config pangraph_groups
# groups; each group is built as a subgraph, and the subgraphs are joined two at
# a time by `pangraph merge` in a balanced tree. MG1655 is only in group 0, so
# the final graph holds one copy of it; group 0 stays on the left (the base) of
# every merge.
PANGRAPH_GROUPS = int(config.get("pangraph_groups", 16))
# Keep in sync with the profile's pangraph_subgraph entry.
PANGRAPH_SUBGRAPH_MEM_MB = 64000
PANGRAPH_SUBGRAPH_REQUEST_MIB = 61035
PANGRAPH_MERGE_MEM_MB_PER_LEVEL = 64000


def pangraph_graph(lo, hi):
    """The graph of groups lo..hi: a group's subgraph, or the merge of its two halves."""
    lo, hi = int(lo), int(hi)
    if lo == hi:
        return PANGRAPH_OUT_DIR / "subgraphs" / f"{lo}.json.gz"
    return PANGRAPH_OUT_DIR / "merges" / f"{lo}-{hi}.json.gz"


def pangraph_halves(wildcards):
    lo, hi = int(wildcards.lo), int(wildcards.hi)
    mid = (lo + hi) // 2
    return pangraph_graph(lo, mid), pangraph_graph(mid + 1, hi)


def pangraph_merge_mem_mb(wildcards):
    """Memory for a merge, growing with the log2 of the number of groups it joins."""
    span = int(wildcards.hi) - int(wildcards.lo) + 1
    return min(PANGRAPH_MERGE_MEM_MB_PER_LEVEL * math.ceil(math.log2(span)), 512000)


PANGRAPH_GRAPH = pangraph_graph(0, PANGRAPH_GROUPS - 1)

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
        fasta = expand(PANGRAPH_OUT_DIR / "input" / "{group}.fa.gz", group=range(PANGRAPH_GROUPS)),
        samples = PANGRAPH_OUT_DIR / "input_samples.txt",
    log: PANGRAPH_LOGS_DIR / "pangraph_input.log"
    benchmark: BENCHMARKS_DIR / "pangraph_input.tsv"
    params:
        ref_name = lookup(dpath="pangraph_ref_name", within=config, default="MG1655"),
        groups = PANGRAPH_GROUPS,
        input_dir = str(PANGRAPH_OUT_DIR / "input"),
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        # The script sorts the assemblies by name (PanGraph's graph depends on
        # input order) and splits them into contiguous groups.
        r"""
        printf '%s\n' {input.samples:q} > {output.samples:q}
        python {PREPARE_PANGRAPH_INPUT:q} --reference {input.reference:q} --ref-name {params.ref_name:q} \
            --groups {params.groups} --out {params.input_dir:q}/'{{group}}.fa.gz' @{output.samples:q} > {log:q} 2>&1
        """

rule pangraph_subgraph:
    input: PANGRAPH_OUT_DIR / "input" / "{group}.fa.gz"
    output: PANGRAPH_OUT_DIR / "subgraphs" / "{group}.json.gz"
    log: PANGRAPH_LOGS_DIR / "pangraph_subgraph" / "{group}.log"
    benchmark: BENCHMARKS_DIR / "pangraph_subgraph_{group}.tsv"
    wildcard_constraints:
        group = r"\d+",
    container: CONTAINERS.format(PANGRAPH_IMAGE)
    # Locally, 10 threads kept only about 2.8 cores busy and 4 were 7% slower.
    threads: 8
    resources:
        mem_mb = PANGRAPH_SUBGRAPH_MEM_MB,
        htcondor_request_mem_mb = PANGRAPH_SUBGRAPH_REQUEST_MIB,
    shell:
        # --verify checks that every input sequence is rebuilt exactly from the graph.
        r"""
        pangraph --verbosity info build {input:q} -o {output:q} -j {threads} --verify --no-progress-bar > {log:q} 2>&1
        """

rule pangraph_merge:
    input:
        left = lambda wc: pangraph_halves(wc)[0],
        right = lambda wc: pangraph_halves(wc)[1],
    output: PANGRAPH_OUT_DIR / "merges" / "{lo}-{hi}.json.gz"
    log: PANGRAPH_LOGS_DIR / "pangraph_merge" / "{lo}-{hi}.log"
    benchmark: BENCHMARKS_DIR / "pangraph_merge_{lo}-{hi}.tsv"
    wildcard_constraints:
        lo = r"\d+",
        hi = r"\d+",
    container: CONTAINERS.format(PANGRAPH_IMAGE)
    threads: 8
    resources:
        mem_mb = pangraph_merge_mem_mb,
        htcondor_request_mem_mb = lambda wc: pangraph_merge_mem_mb(wc) * 1000 * 1000 // 1024 ** 2,
    shell:
        r"""
        pangraph --verbosity info merge {input.left:q} {input.right:q} -o {output:q} -j {threads} --verify > {log:q} 2>&1
        """

rule pangraph_features:
    input: PANGRAPH_GRAPH
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

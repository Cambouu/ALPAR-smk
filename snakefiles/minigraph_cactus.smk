from pathlib import Path

# Minigraph-Cactus experiment: variants from a pangenome graph of the assemblies themselves,
# written by vg deconstruct as a VCF in reference (MG1655) coordinates, in place of snippy's
# SNPs. Minigraph-Cactus 3.3.0 builds one graph of every strain plus the reference, without
# labels; cactus-pangenome's --vcf step runs vg deconstruct on the MG1655 path, so each strain's
# genotype is the path its own assembly takes through each bubble (no reads, no pseudo-reads).
# mc_vcf_to_features.py splits and atomises the VCF with bcftools and writes ALPAR's long
# feature table (mc|snp|<pos>|<REF>><ALT>, mc|ind|...).
#
# Each set in config mc_feature_sets (default GPA + MC) is merged, pivoted once for all
# antibiotics and trained with the baseline combined_ml. The SNP + GPA and GPA + PanGraph
# results of experiment/pangraph (same split, same code) are reported as reference rows and are
# not retrained. Not part of automatix; target experiment_minigraph_cactus. Every rule here
# writes under minigraph_cactus/ or experiments/minigraph_cactus/.

MC_OUT_DIR = OUT_DIR / "minigraph_cactus"
MC_LOGS_DIR = MC_OUT_DIR / "logs"
EXPERIMENT_MC_DIR = OUT_DIR / "experiments" / "minigraph_cactus"

PREPARE_MC_INPUT = SCRIPTS_DIR / "prepare_mc_input.py"
MC_VCF_TO_FEATURES = SCRIPTS_DIR / "mc_vcf_to_features.py"

# alpar-smk-minigraph-cactus:1.0.0: the official cactus v3.3.0 image (cactus-pangenome, vg,
# bcftools) plus the Snakemake worker (minigraph-cactus/containers/minigraph-cactus.Dockerfile in
# the thesis repository). The bioconda cactus package is incomplete, so the rules have no conda env.
MC_IMAGE = "minigraph-cactus@sha256:2116d2cf79c781c89040ca403356d04e1385fcd4e8f84f7c8342477b563560aa"

# HTCondor's docker universe does not pass on the image's ENV; the Cactus tools need these.
MC_ENV = (
    "export PATH=/usr/local/cactus/cactus_env/bin:/usr/local/cactus/bin:$PATH"
    " PYTHONPATH=/usr/local/cactus LD_LIBRARY_PATH=/usr/local/cactus/lib CACTUS_INSIDE_CONTAINER=1"
)

MC_FEATURE_SETS = config.get("mc_feature_sets", {
    "gpa_mc": ["gpa", "mc"],
})

wildcard_constraints:
    mc_set="(?:" + "|".join(re.escape(name) for name in MC_FEATURE_SETS) + ")"

# Same model settings as the baseline target ml.
MC_ML_SETTINGS = dict(
    random_seed=[42],
    test_size=[0.2],
    resampling_strategy=["cv"],
    model_type=["xgb"],
    feature_importance_analysis_strategy=["gini"],
)

# Keep in sync with the profile's mc_build entry.
MC_BUILD_THREADS = int(config.get("mc_build_threads", 64))
MC_BUILD_MEM_MB = int(config.get("mc_build_mem_mb", 512000))
# Which deconstruct VCF becomes the features: "full" (no clipping) or "clip" (cactus default),
# optionally with the vcfwave-realigned version ("full.wave", "wave" for clip).
MC_VCF = config.get("mc_vcf", "full")


def mc_vcf_path(kind):
    stem = {"clip": "mc", "wave": "mc.wave"}.get(kind, f"mc.{kind}")
    return MC_OUT_DIR / "graph" / f"{stem}.vcf.gz"


# -----------------------
# Graph and Features
# -----------------------

rule mc_input:
    input:
        store = rules.rename_files.output.store,
        # Every strain of the run: the graph is built without labels.
        samples = lambda wc: expand(
            Path(rules.rename_files.output.store) / "{sample}",
            sample = get_sample_names(wc)
        ),
        # MG1655 as GenBank, the reference snippy uses, so positions share its coordinates.
        reference = GBFF_FILE,
    output:
        seqfile = MC_OUT_DIR / "input" / "seqfile.txt",
        reference = MC_OUT_DIR / "input" / "MG1655.fa",
        samples = MC_OUT_DIR / "input_samples.txt",
    log: MC_LOGS_DIR / "mc_input.log"
    benchmark: BENCHMARKS_DIR / "mc_input.tsv"
    params:
        ref_name = lookup(dpath="mc_ref_name", within=config, default="MG1655"),
        input_dir = str(MC_OUT_DIR / "input"),
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        # Uppercased copies with short contig names; the reference first, then the assemblies sorted.
        r"""
        printf '%s\n' {input.samples:q} > {output.samples:q}
        python {PREPARE_MC_INPUT:q} --reference {input.reference:q} --ref-name {params.ref_name:q} \
            --outdir {params.input_dir:q} --seqfile {output.seqfile:q} @{output.samples:q} > {log:q} 2>&1
        """

rule mc_build:
    input:
        seqfile = rules.mc_input.output.seqfile,
        reference = rules.mc_input.output.reference,
    output:
        full = mc_vcf_path("full"),
        clip = mc_vcf_path("clip"),
        full_wave = mc_vcf_path("full.wave"),
        clip_wave = mc_vcf_path("wave"),
        gbz = MC_OUT_DIR / "graph" / "mc.full.gbz",
        pathlen = MC_OUT_DIR / "graph" / "pathlen.full.tsv",
    log: MC_LOGS_DIR / "mc_build.log"
    benchmark: BENCHMARKS_DIR / "mc_build.tsv"
    params:
        ref_name = lookup(dpath="mc_ref_name", within=config, default="MG1655"),
        out_dir = str(MC_OUT_DIR / "graph"),
        # --noSplit: one graph for the whole genome. The default splits the assemblies by reference
        # chromosome, and with MG1655's single chromosome it dropped about twice as much sequence
        # (plasmids and draft contigs that fail the coverage filters) in the local 2-genome check.
        extra = lookup(dpath="mc_extra_args", within=config, default="--noSplit"),
        max_mem = lambda wc, resources: f"{int(resources.mem_mb * 0.94) // 1000}G",
    container: CONTAINERS.format(MC_IMAGE)
    threads: MC_BUILD_THREADS
    resources:
        mem_mb = MC_BUILD_MEM_MB,
        htcondor_request_mem_mb = MC_BUILD_MEM_MB * 1000 * 1000 // 1024 ** 2,
    shell:
        # minGAFBlockLength=0 keeps draft contigs >= 250 kb whose longest minigraph block is shorter
        # (docs/vg-pilot-methods.md in the thesis repository). The Toil job store and work files
        # stay on the job's node-local /tmp; cactus writes the outputs to out_dir.
        r"""
        {MC_ENV}
        cfg=$(python3 -c 'import cactus,os;print(os.path.join(os.path.dirname(cactus.__file__),"cactus_progressive_config.xml"))')
        sed 's/minGAFBlockLength="250000"/minGAFBlockLength="0"/' "$cfg" > "$TMPDIR/config.xml"
        grep -q 'minGAFBlockLength="0"' "$TMPDIR/config.xml"
        rm -rf "$TMPDIR/js"
        cactus-pangenome "$TMPDIR/js" {input.seqfile:q} --outDir {params.out_dir:q} --outName mc \
            --reference {params.ref_name:q} --configFile "$TMPDIR/config.xml" \
            --vcf clip full --vcfwave --gbz clip full --gfa clip full {params.extra} \
            --maxCores {threads} --maxMemory {params.max_mem} --workDir "$TMPDIR" > {log:q} 2>&1
        vg paths -E -x {output.gbz:q} > {output.pathlen:q} 2>> {log:q}
        """

rule mc_features:
    input:
        vcf = mc_vcf_path(MC_VCF),
        reference = rules.mc_input.output.reference,
    output: MC_OUT_DIR / "mc_features.tsv"
    log: MC_LOGS_DIR / "mc_features.log"
    benchmark: BENCHMARKS_DIR / "mc_features.tsv"
    params:
        ref_name = lookup(dpath="mc_ref_name", within=config, default="MG1655"),
    container: CONTAINERS.format(MC_IMAGE)
    threads: 1
    resources:
        mem_mb = 32000,
        htcondor_request_mem_mb = 30518,
    shell:
        r"""
        {MC_ENV}
        python3 {MC_VCF_TO_FEATURES:q} {input.vcf:q} {input.reference:q} {output:q} \
            --ref-sample {params.ref_name:q} > {log:q} 2>&1
        """

rule minigraph_cactus:
    localrule: True
    input: rules.mc_features.output
    output: touch(OUT_DIR / "flags" / "minigraph_cactus.done")


# -----------------------
# Feature Table and ML
# -----------------------

MC_FAMILY_TABLES = {
    "snp": rules.binary_mutation_table.output[0],
    "gpa": rules.binary_gpa.output[0],
    "pangraph": rules.pangraph_features.output.features,
    "mc": rules.mc_features.output[0],
}

# All families are built from every strain, so one table serves every antibiotic;
# combined_ml takes each antibiotic's strains from its split.
rule merge_mc_set:
    localrule: True
    input: lambda wc: [MC_FAMILY_TABLES[family] for family in MC_FEATURE_SETS[wc.mc_set]]
    output: temp(EXPERIMENT_MC_DIR / "{mc_set}" / "feature_table" / "merged_table.tsv"),
    threads: 1,
    shell:
        r"""
        cat {input:q} > {output:q}
        """

use rule pivot_merged_features_miller as pivot_mc_set with:
    group: "pivot_mc_set_batch"
    input: rules.merge_mc_set.output,
    output: EXPERIMENT_MC_DIR / "{mc_set}" / "feature_table" / "merged_table_pivot.tsv",
    log: EXPERIMENT_MC_DIR / "{mc_set}" / "feature_table" / "logs" / "pivot_mc_set.log",
    benchmark: BENCHMARKS_DIR / "pivot_mc_set_{mc_set}.tsv",


MC_ML_STEM = str(EXPERIMENT_MC_DIR / "{mc_set}" / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}")

use rule combined_ml as combined_ml_mc with:
    group: "combined_ml_mc_batch"
    input:
        binary_mutation_table = rules.pivot_mc_set.output[0],
        phenotype_table = rules.phenotype_dataframe_creator.output[0],
        train = lambda wildcards: expand(rules.split_train_test.output, split_category=["train"], antibiotic=wildcards.antibiotic)[0],
        test = lambda wildcards: expand(rules.split_train_test.output, split_category=["test"], antibiotic=wildcards.antibiotic)[0],
    output:
        best_params = MC_ML_STEM + "_best_params.txt",
        model_file  = MC_ML_STEM + "_model.sav",
        result      = MC_ML_STEM + "_Result.txt",
        fia         = MC_ML_STEM + ".txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_mc_{mc_set}_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.tsv"
    log: MC_ML_STEM + ".log"
    params:
        feature_importance_analysis = True,
        save_model = True,
        param_grid_size = lookup(dpath="ml_param_grid_size", within=config, default="small"),
        # As in the PanGraph experiment, so the reference rows are comparable.
        param_grid_low_memory_mode = True,
        parallel_fits = lookup(dpath="ml_parallel_fits", within=config, default=1),


# -----------------------
# Report
# -----------------------

def mc_ml_outputs(name):
    return expand(
        getattr(rules.combined_ml_mc.output, name),
        mc_set=list(MC_FEATURE_SETS),
        antibiotic=ANTIBIOTICS,
        **MC_ML_SETTINGS,
    )

rule mc_report:
    input:
        results = mc_ml_outputs("result"),
        best_params = mc_ml_outputs("best_params"),
        fia = mc_ml_outputs("fia"),
        pivots = expand(rules.pivot_mc_set.output[0], mc_set=list(MC_FEATURE_SETS)),
        features = rules.mc_features.output[0],
        samples = rules.mc_input.output.samples,
        pathlen = rules.mc_build.output.pathlen,
        seqfile = rules.mc_input.output.seqfile,
    output:
        markdown = EXPERIMENT_MC_DIR / "report.md",
        table = EXPERIMENT_MC_DIR / "report.tsv",
    log: MC_LOGS_DIR / "mc_report.log"
    params:
        feature_sets = MC_FEATURE_SETS,
        vcf = MC_VCF,
        features_log = rules.mc_features.log[0],
        # Earlier runs on the same split, shown for reference; not rerun.
        reference_ml_dirs = {
            "SNP + GPA (experiment/pangraph)": str(OUT_DIR / "experiments" / "pangraph" / "snp_gpa" / "ml"),
            "GPA + PanGraph (experiment/pangraph)": str(OUT_DIR / "experiments" / "pangraph" / "gpa_pangraph" / "ml"),
        },
        benchmarks_dir = str(BENCHMARKS_DIR),
        # Quinolone resistance genes, MG1655 (NC_000913.3) coordinates, both on the minus strand.
        loci = {"gyrA": [2336793, 2339420], "parC": [3163715, 3165973]},
        top_features = 20,
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    resources:
        mem_mb = 16000,
    script:
        SCRIPTS_DIR / "mc_report.py"

rule experiment_minigraph_cactus:
    localrule: True
    input: rules.mc_report.output
    output: touch(OUT_DIR / "flags" / "experiment_minigraph_cactus.done")

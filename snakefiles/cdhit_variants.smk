import re
from pathlib import Path

# CD-HIT variant experiment: do clusters that merge more proteins change the
# models? Each variant reruns the GPA step's CD-HIT clustering (rule
# cdhit_runner) with some options changed, and its clusters replace the
# baseline's in both the GPA table and the panproteome: protein graphs,
# BubbleGun bubbles, bubble features (the fixed command of panproteome_fixed.smk)
# and the PanPA test alignments. Two feature sets per variant are trained with
# the baseline combined_ml:
#   snp_gpa              SNP + GPA; compare experiments/snp_gpa
#   snp_gpa_panproteome  SNP + GPA + panproteome; compare
#                        experiments/snp_gpa_panproteome_fixed
# The SNP table, phenotypes and DataSAIL splits are the baseline's. Not part of
# automatix; target experiment_cdhit_variants. Every rule here writes under
# cdhit_variants/ or experiments/cdhit_variants/.

CDHIT_VARIANTS_OUT_DIR = OUT_DIR / "cdhit_variants"
CDHIT_VARIANTS_LOGS_DIR = CDHIT_VARIANTS_OUT_DIR / "logs"
EXPERIMENT_CDHIT_DIR = OUT_DIR / "experiments" / "cdhit_variants"
MERGE_CDHIT_STEPS = SCRIPTS_DIR / "merge_cdhit_steps.py"

# Options of rule cdhit_runner, plus cd-hit's defaults for the ones it leaves
# unset, so a variant names only what it changes.
CDHIT_RUNNER_OPTIONS = {
    "c": 0.7, "n": 5, "G": 1, "s": 0.0, "aL": 0.0, "AL": 99999999,
    "aS": 0.0, "AS": 99999999, "l": 10, "b": 20, "M": 0, "d": 0,
}

# Cluster counts from the CD-HIT sweep of 2026-11-07 (HTCondor 71972; thesis
# repository outputs/gene-frequency-2026-11-07/cdhit-sweep), over the 1,814
# genomes with at least 1,000 proteins. Baseline: 34,062 clusters, 10,106 in a
# single genome.
CDHIT_VARIANTS = config.get("cdhit_variants", {
    # Local identity (identical residues over the alignment length) on at least
    # 70% of the shorter protein; proteins under 50 aa are left out (cd-hit -l).
    # 31,777 clusters, 9,165 in a single genome.
    "G0_aS0.7_l50": {"options": {"G": 0, "aS": 0.7, "l": 50}, "two_step": False},
    # The baseline representatives clustered again at 0.6 identity (two-step
    # clustering from the CD-HIT guide). 29,131 clusters, 8,009 in a single genome.
    "hier_c0.6_n4": {"options": {"c": 0.6, "n": 4}, "two_step": True},
})
CDHIT_VARIANT_SETS = ("snp_gpa", "snp_gpa_panproteome")
CDHIT_TWO_STEP_VARIANTS = [name for name, spec in CDHIT_VARIANTS.items() if spec.get("two_step")]

wildcard_constraints:
    cdhit_variant = "(?:" + "|".join(re.escape(name) for name in CDHIT_VARIANTS) + ")",
    cdhit_set = "(?:" + "|".join(CDHIT_VARIANT_SETS) + ")",

CDHIT_VARIANT_DIR = CDHIT_VARIANTS_OUT_DIR / "{cdhit_variant}"
CDHIT_VARIANT_LOGS = CDHIT_VARIANTS_LOGS_DIR / "{cdhit_variant}"
CDHIT_VARIANT_PANGENOME = CDHIT_VARIANT_DIR / "pangenome"
CDHIT_VARIANT_SHARD_MANIFEST = CDHIT_VARIANT_PANGENOME / "cluster_shards" / "{shard}.txt"


def is_two_step(wildcards):
    return bool(CDHIT_VARIANTS[wildcards.cdhit_variant].get("two_step"))


def cdhit_variant_options(wildcards):
    options = {**CDHIT_RUNNER_OPTIONS, **CDHIT_VARIANTS[wildcards.cdhit_variant]["options"]}
    return " ".join(f"-{key} {value}" for key, value in options.items())


# -----------------------
# Clusters
# -----------------------

rule cdhit_variant:
    group: "cdhit_variant_batch"
    input:
        # A two-step variant clusters the baseline representatives again.
        lambda wc: rules.cdhit_runner.output.faa if is_two_step(wc) else rules.combine_faa_files.output[0],
    output:
        faa = CDHIT_VARIANT_DIR / "cd-hit" / "cdhit_output.faa",
        clstr = CDHIT_VARIANT_DIR / "cd-hit" / "cdhit_output.faa.clstr",
    log: CDHIT_VARIANT_LOGS / "cdhit_variant.log"
    benchmark: BENCHMARKS_DIR / "cdhit_variant_{cdhit_variant}.tsv"
    params:
        options = cdhit_variant_options,
    conda: ENVS_DIR.format("cd-hit")
    container: CONTAINERS.format("cd-hit:1.0.0")
    threads: 8
    shell:
        r"""
        cd-hit -i {input:q} -o {output.faa:q} -T {threads} {params.options} > {log:q} 2>&1
        """


rule merge_cdhit_steps:
    group: "merge_cdhit_steps_batch"
    input:
        first = rules.cdhit_runner.output.clstr,
        second = rules.cdhit_variant.output.clstr,
    output: CDHIT_VARIANT_DIR / "cd-hit" / "two_step.clstr"
    log: CDHIT_VARIANT_LOGS / "merge_cdhit_steps.log"
    benchmark: BENCHMARKS_DIR / "merge_cdhit_steps_{cdhit_variant}.tsv"
    wildcard_constraints:
        cdhit_variant = "(?:" + "|".join(re.escape(name) for name in CDHIT_TWO_STEP_VARIANTS) + ")",
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    shell:
        r"""
        python {MERGE_CDHIT_STEPS:q} {input.first:q} {input.second:q} {output:q} 2> {log:q}
        """


def cdhit_variant_clusters(wildcards):
    """The variant's final .clstr: the merged file for two-step variants."""
    rule_output = rules.merge_cdhit_steps.output[0] if is_two_step(wildcards) else rules.cdhit_variant.output.clstr
    return expand(rule_output, cdhit_variant=wildcards.cdhit_variant)[0]


use rule binary_gpa_cdhit as binary_gpa_cdhit_variant with:
    group: "binary_gpa_cdhit_variant_batch"
    input: cdhit_variant_clusters,
    output: CDHIT_VARIANT_DIR / "binary_gpa_cdhit.tsv"
    benchmark: BENCHMARKS_DIR / "binary_gpa_cdhit_variant_{cdhit_variant}.tsv"
    log: CDHIT_VARIANT_LOGS / "binary_gpa_cdhit.log"


# -----------------------
# Panproteome
# -----------------------

use rule split_cluster_fasta as split_cluster_fasta_variant with:
    group: "split_cluster_fasta_variant_batch"
    input:
        cdhit_clstr = cdhit_variant_clusters,
        combined_proteins = rules.combine_faa_files.output[0],
    output: directory(CDHIT_VARIANT_PANGENOME / "cluster_sequences"),
    log: CDHIT_VARIANT_LOGS / "split_cluster_fasta.log",
    benchmark: BENCHMARKS_DIR / "split_cluster_fasta_variant_{cdhit_variant}.tsv",


use rule shard_clusters as shard_clusters_variant with:
    input: rules.split_cluster_fasta_variant.output[0],
    output: expand(CDHIT_VARIANT_SHARD_MANIFEST, shard=PANGENOME_SHARDS, allow_missing=True),


use rule cluster_fasta_splits as cluster_fasta_splits_variant with:
    group: "cluster_fasta_splits_variant_batch"
    input:
        cluster_store = rules.split_cluster_fasta_variant.output[0],
        datasail_splits = rules.datasail_runner.output[0],
    output: directory(CDHIT_VARIANT_PANGENOME / "cluster_fasta_splits" / "{antibiotic}" / "{split_category}"),
    log: CDHIT_VARIANT_LOGS / "cluster_fasta_splits" / "{antibiotic}_{split_category}.log",
    benchmark: BENCHMARKS_DIR / "cluster_fasta_splits_variant_{cdhit_variant}_{antibiotic}_{split_category}.tsv",


use rule align_clusters as align_clusters_variant with:
    input:
        manifest = CDHIT_VARIANT_SHARD_MANIFEST,
        cluster_store = rules.split_cluster_fasta_variant.output[0],
    output: directory(CDHIT_VARIANT_PANGENOME / "cluster_alignments" / "{shard}")
    log: CDHIT_VARIANT_LOGS / "align_clusters" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "align_clusters_variant_{cdhit_variant}_{shard}.tsv"


use rule panpa_build_gfa as panpa_build_gfa_variant with:
    input:
        manifest = CDHIT_VARIANT_SHARD_MANIFEST,
        alignments = rules.align_clusters_variant.output[0],
    output: directory(CDHIT_VARIANT_PANGENOME / "panpa" / "gfa" / "{shard}")
    log: CDHIT_VARIANT_LOGS / "panpa_build_gfa" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "panpa_build_gfa_variant_{cdhit_variant}_{shard}.tsv"


use rule bubblegun_runner as bubblegun_runner_variant with:
    input:
        manifest = CDHIT_VARIANT_SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa_variant.output[0],
    output: directory(CDHIT_VARIANT_PANGENOME / "bubblegun" / "{shard}")
    log: CDHIT_VARIANT_LOGS / "bubblegun_runner" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubblegun_runner_variant_{cdhit_variant}_{shard}.tsv"


# The fixed bubble_features command (path cap, warning-level logs), as in
# panproteome_fixed.smk, so the panproteome sets compare like for like.
use rule bubble_features_nt as bubble_features_variant with:
    input:
        manifest = CDHIT_VARIANT_SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa_variant.output[0],
        bubblegun_dir = rules.bubblegun_runner_variant.output[0],
        phenotype_table = lambda wc: rules.split_phenotype_dataframe.output[0].format(
            antibiotic=wc.antibiotic, split_category="train"
        ),
    output:
        output_dir = directory(CDHIT_VARIANT_PANGENOME / "bubble_features" / "train" / "{antibiotic}" / "{shard}"),
        lor_lookup_dir = directory(CDHIT_VARIANT_PANGENOME / "bubble_features_{antibiotic}_lor_lookup" / "{shard}"),
    log: CDHIT_VARIANT_LOGS / "bubble_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubble_features_variant_{cdhit_variant}_{antibiotic}_{shard}.tsv"


# Merges fail on an empty table, like the nucleotide and fixed ones.
use rule merge_bubble_features_nt as merge_bubble_features_variant with:
    input: expand(rules.bubble_features_variant.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: CDHIT_VARIANT_PANGENOME / "bubble_features_{antibiotic}.tsv"


use rule panpa_align as panpa_align_variant with:
    input:
        manifest = CDHIT_VARIANT_SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa_variant.output[0],
        test_sequence_folder = lambda wc: expand(
            rules.cluster_fasta_splits_variant.output[0],
            cdhit_variant=wc.cdhit_variant, antibiotic=wc.antibiotic, split_category="test",
        )[0],
    output: directory(CDHIT_VARIANT_PANGENOME / "panpa" / "alignments" / "{antibiotic}" / "{shard}")
    log: CDHIT_VARIANT_LOGS / "panpa_align" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "panpa_align_variant_{cdhit_variant}_{antibiotic}_{shard}.tsv"


use rule gaf_lor_features as gaf_lor_features_variant with:
    input:
        manifest = CDHIT_VARIANT_SHARD_MANIFEST,
        gaf_dir = rules.panpa_align_variant.output[0],
        lor_lookup_dir = rules.bubble_features_variant.output.lor_lookup_dir,
    output:
        output_dir = directory(CDHIT_VARIANT_PANGENOME / "bubble_features" / "test" / "{antibiotic}" / "{shard}"),
    log: CDHIT_VARIANT_LOGS / "gaf_lor_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "gaf_lor_features_variant_{cdhit_variant}_{antibiotic}_{shard}.tsv"


use rule merge_bubble_features_nt as merge_bubble_features_test_variant with:
    input: expand(rules.gaf_lor_features_variant.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: CDHIT_VARIANT_PANGENOME / "bubble_features_test_{antibiotic}.tsv"


# -----------------------
# Feature Tables and ML
# -----------------------

# SNP + GPA holds every strain, so one table serves both antibiotics; combined_ml
# takes each antibiotic's strains from its split.
rule merge_cdhit_variant_snp_gpa:
    localrule: True
    input:
        snp = rules.binary_mutation_table.output[0],
        gpa = rules.binary_gpa_cdhit_variant.output[0],
    output: temp(EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "snp_gpa" / "feature_table" / "merged_table.tsv"),
    shell:
        r"""
        cat {input.snp:q} {input.gpa:q} > {output:q}
        """


# One pivot per job: a bundled SNP + GPA pivot would wait for the panproteome.
use rule pivot_merged_features_miller as pivot_cdhit_variant_snp_gpa with:
    group: "pivot_cdhit_variant_snp_gpa_batch"
    input: rules.merge_cdhit_variant_snp_gpa.output,
    output: EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "snp_gpa" / "feature_table" / "merged_table_pivot.tsv",
    log: EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "snp_gpa" / "feature_table" / "logs" / "pivot.log",
    benchmark: BENCHMARKS_DIR / "pivot_cdhit_variant_snp_gpa_{cdhit_variant}.tsv",


use rule merge_features as merge_features_cdhit_variant with:
    input:
        snp = rules.binary_mutation_table.output,
        gpa = rules.binary_gpa_cdhit_variant.output,
        pangenome_train = rules.merge_bubble_features_variant.output,
        pangenome_test = rules.merge_bubble_features_test_variant.output,
    output: temp(EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "snp_gpa_panproteome" / "feature_table" / "merged_table_{antibiotic}.tsv"),


use rule pivot_merged_features_miller as pivot_cdhit_variant_panproteome with:
    group: "pivot_cdhit_variant_panproteome_batch"
    input: rules.merge_features_cdhit_variant.output,
    output: EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "snp_gpa_panproteome" / "feature_table" / "merged_table_pivot_{antibiotic}.tsv",
    log: EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "snp_gpa_panproteome" / "feature_table" / "logs" / "pivot_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_cdhit_variant_panproteome_{cdhit_variant}_{antibiotic}.tsv",


def cdhit_variant_feature_table(wildcards):
    if wildcards.cdhit_set == "snp_gpa":
        pattern = rules.pivot_cdhit_variant_snp_gpa.output[0]
    else:
        pattern = rules.pivot_cdhit_variant_panproteome.output[0]
    return expand(pattern, cdhit_variant=wildcards.cdhit_variant, antibiotic=wildcards.antibiotic)[0]


# Same model, grid and search as the baseline combined_ml.
CDHIT_VARIANT_ML_STEM = str(EXPERIMENT_CDHIT_DIR / "{cdhit_variant}" / "{cdhit_set}" / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}")

use rule combined_ml as combined_ml_cdhit_variant with:
    group: "combined_ml_cdhit_variant_batch"
    input:
        binary_mutation_table = cdhit_variant_feature_table,
        phenotype_table = rules.phenotype_dataframe_creator.output[0],
        train = lambda wildcards: expand(rules.split_train_test.output, split_category=["train"], antibiotic=wildcards.antibiotic)[0],
        test = lambda wildcards: expand(rules.split_train_test.output, split_category=["test"], antibiotic=wildcards.antibiotic)[0],
    output:
        best_params = CDHIT_VARIANT_ML_STEM + "_best_params.txt",
        model_file  = CDHIT_VARIANT_ML_STEM + "_model.sav",
        result      = CDHIT_VARIANT_ML_STEM + "_Result.txt",
        fia         = CDHIT_VARIANT_ML_STEM + ".txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_cdhit_variant_{cdhit_variant}_{cdhit_set}_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.tsv"
    log: CDHIT_VARIANT_ML_STEM + ".log"


rule experiment_cdhit_variants:
    localrule: True
    input:
        expand(
            rules.combined_ml_cdhit_variant.output,
            cdhit_variant=list(CDHIT_VARIANTS),
            cdhit_set=CDHIT_VARIANT_SETS,
            antibiotic=ANTIBIOTICS,
            random_seed=[42],
            test_size=[0.2],
            resampling_strategy=["cv"],
            model_type=["xgb"],
            feature_importance_analysis_strategy=["gini"],
        )
    output: touch(OUT_DIR / "flags" / "experiment_cdhit_variants.done")

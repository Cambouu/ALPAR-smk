from pathlib import Path

# Fixed panproteome baseline: SNP + GPA + panproteome, with bubble features
# recomputed by the corrected bubble_features.py. In the original baseline,
# about 3,800 clusters per antibiotic hit an UnboundLocalError that was
# swallowed after their outputs were opened, so their bubble features were
# truncated (run 67228 investigation, 2026-10-17).
#
# Reuses the baseline protein graphs, BubbleGun output and test alignments,
# which the bug did not affect, and uses the same bubble_features command as
# the nucleotide experiment (path cap and log level), so both compare like
# for like. The original baseline outputs are untouched; everything here
# writes under panproteome_fixed/ or experiments/snp_gpa_panproteome_fixed/.

PANPROTEOME_FIXED_OUT_DIR = OUT_DIR / "panproteome_fixed"
PANPROTEOME_FIXED_LOGS_DIR = PANPROTEOME_FIXED_OUT_DIR / "logs"
EXPERIMENT_FIXED_DIR = OUT_DIR / "experiments" / "snp_gpa_panproteome_fixed"


use rule bubble_features_nt as bubble_features_fixed with:
    input:
        manifest = SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa.output[0],
        bubblegun_dir = rules.bubblegun_runner.output[0],
        phenotype_table = lambda wc: rules.split_phenotype_dataframe.output[0].format(
            antibiotic=wc.antibiotic, split_category="train"
        ),
    output:
        output_dir = directory(PANPROTEOME_FIXED_OUT_DIR / "bubble_features" / "train" / "{antibiotic}" / "{shard}"),
        lor_lookup_dir = directory(PANPROTEOME_FIXED_OUT_DIR / "bubble_features_{antibiotic}_lor_lookup" / "{shard}"),
    log: PANPROTEOME_FIXED_LOGS_DIR / "bubble_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubble_features_fixed_{antibiotic}_{shard}.tsv"


use rule merge_bubble_features as merge_bubble_features_fixed with:
    input: expand(rules.bubble_features_fixed.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: PANPROTEOME_FIXED_OUT_DIR / "bubble_features_{antibiotic}.tsv"


use rule gaf_lor_features as gaf_lor_features_fixed with:
    input:
        manifest = SHARD_MANIFEST,
        gaf_dir = rules.panpa_align.output[0],
        lor_lookup_dir = rules.bubble_features_fixed.output.lor_lookup_dir,
    output:
        output_dir = directory(PANPROTEOME_FIXED_OUT_DIR / "bubble_features" / "test" / "{antibiotic}" / "{shard}"),
    log: PANPROTEOME_FIXED_LOGS_DIR / "gaf_lor_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "gaf_lor_features_fixed_{antibiotic}_{shard}.tsv"


use rule merge_bubble_features_test as merge_bubble_features_test_fixed with:
    input: expand(rules.gaf_lor_features_fixed.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: PANPROTEOME_FIXED_OUT_DIR / "bubble_features_test_{antibiotic}.tsv"


# -----------------------
# Feature Table and ML
# -----------------------

use rule merge_features as merge_features_fixed with:
    input:
        snp = rules.binary_mutation_table.output,
        gpa = rules.binary_gpa.output,
        pangenome_train = rules.merge_bubble_features_fixed.output,
        pangenome_test = rules.merge_bubble_features_test_fixed.output,
    output: EXPERIMENT_FIXED_DIR / "feature_table" / "merged_table_{antibiotic}.tsv",


use rule pivot_merged_features_miller as pivot_features_fixed with:
    group: "pivot_features_fixed_batch"
    input: rules.merge_features_fixed.output,
    output: EXPERIMENT_FIXED_DIR / "feature_table" / "merged_table_pivot_{antibiotic}.tsv",
    log: EXPERIMENT_FIXED_DIR / "feature_table" / "logs" / "pivot_features_fixed_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_features_fixed_{antibiotic}.tsv",


EXPERIMENT_FIXED_ML_STEM = str(EXPERIMENT_FIXED_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}")

use rule combined_ml as combined_ml_fixed with:
    group: "combined_ml_fixed_batch"
    input:
        binary_mutation_table = rules.pivot_features_fixed.output[0],
        phenotype_table = rules.phenotype_dataframe_creator.output[0],
        train = lambda wildcards: expand(rules.split_train_test.output, split_category=["train"], **wildcards)[0],
        test = lambda wildcards: expand(rules.split_train_test.output, split_category=["test"], **wildcards)[0],
    output:
        best_params = EXPERIMENT_FIXED_ML_STEM + "_best_params.txt",
        model_file  = EXPERIMENT_FIXED_ML_STEM + "_model.sav",
        result      = EXPERIMENT_FIXED_ML_STEM + "_Result.txt",
        fia         = EXPERIMENT_FIXED_ML_STEM + ".txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_fixed_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.tsv"
    log: EXPERIMENT_FIXED_ML_STEM + ".log"


rule experiment_panproteome_fixed:
    localrule: True
    input:
        expand(
            rules.combined_ml_fixed.output,
            antibiotic=ANTIBIOTICS,
            random_seed=[42],
            test_size=[0.2],
            resampling_strategy=["cv"],
            model_type=["xgb"],
            feature_importance_analysis_strategy=["gini"]
        )
    output: touch(OUT_DIR / "flags" / "experiment_panproteome_fixed.done")

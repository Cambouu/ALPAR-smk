from pathlib import Path

PRPS_OUT_DIR = OUT_DIR / "prps"
PRPS_LOGS_DIR = PRPS_OUT_DIR / "logs"

rule prps_runner:
    group: "prps_runner_batch"
    input:
        phylogeny_tree = rules.mashtree_runner.output[0],
        feature_matrix = rules.pivot_merged_features_miller.output[0],
    output: PRPS_OUT_DIR / "prps_scores_{antibiotic}.tsv",
    benchmark: BENCHMARKS_DIR / "prps_runner_{antibiotic}.tsv",
    log: PRPS_LOGS_DIR / "prps_runner_{antibiotic}.log",
    conda: ENVS_DIR.format("prps"),
    container: CONTAINERS.format("prps:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "prps.py"

rule prps:
    localrule: True
    input: lambda wc: expand(rules.prps_runner.output, sample = get_sample_names(wc), antibiotic = ANTIBIOTICS)
    output: touch(OUT_DIR / "flags" / "prps.done")


# ALPAR vs ALPAR-smk comparison: PRPS on a feature set's wide table (upstream ALPAR
# computes PRPS on its SNP + GPA table, not on the panproteome-including baseline).
# Same script as prps_runner. The rule-level container overrides the profile default,
# which only names prps_runner. prps.py keeps the whole table as Python floats
# (about 46 GB for 1,868 x 778k), hence the memory request.
use rule prps_runner as prps_runner_feature_set with:
    group: "prps_runner_feature_set_batch"
    input:
        phylogeny_tree = rules.mashtree_runner.output[0],
        feature_matrix = rules.pivot_feature_set.output[0],
    output: EXPERIMENT_DIR / "prps" / "prps_scores_{antibiotic}.tsv",
    benchmark: BENCHMARKS_DIR / "prps_runner_{feature_set}_{antibiotic}.tsv",
    log: EXPERIMENT_DIR / "prps" / "logs" / "prps_runner_{antibiotic}.log",
    resources:
        mem_mb = 96000,
        htcondor_request_mem_mb = 91553,
        container_image = "docker.io/cambouu/alpar-smk-prps:1.0.0",

from pathlib import Path

FEATURE_TABLE_OUT_DIR = OUT_DIR / "feature_table"
FEATURE_TABLES_LOGS_DIR = FEATURE_TABLE_OUT_DIR / "logs"


# ML feature-set experiments. Each entry of config["feature_sets"] enables some
# feature families with use_snp, use_gpa and use_pangenome. Only the merge, pivot
# and ML steps run per set, under experiments/<feature_set>/. They reuse the
# existing feature tables, so no upstream feature computation is repeated.
FEATURE_SETS = config.get("feature_sets", {})
EXPERIMENT_DIR = OUT_DIR / "experiments" / "{feature_set}"

wildcard_constraints:
    feature_set="(?:" + "|".join(re.escape(name) for name in FEATURE_SETS) + ")"


def use_features(family):
    return lookup(dpath=f"feature_sets/{{feature_set}}/use_{family}", within=config, default=False)


rule merge_features:
    localrule: True
    input:
        branch(use_features("snp"), then=rules.binary_mutation_table.output),
        branch(use_features("gpa"), then=rules.binary_gpa.output),
        branch(use_features("pangenome"), then=rules.merge_bubble_features.output),
        branch(use_features("pangenome"), then=rules.merge_bubble_features_test.output),
    output: EXPERIMENT_DIR / "feature_table" / "merged_table_{antibiotic}.tsv",
    threads: 1,
    shell:
        r"""
        cat {input} > {output}
        """

# Baseline: all feature families, at the original paths used by PRPS and ML.
use rule merge_features as merge_all_features with:
    input:
        rules.binary_mutation_table.output,
        rules.binary_gpa.output,
        rules.merge_bubble_features.output,
        rules.merge_bubble_features_test.output,
    output: FEATURE_TABLE_OUT_DIR / "merged_table_{antibiotic}.tsv",

rule pivot_merged_features_miller:
    group: "pivot_merged_features_miller_batch"
    input: rules.merge_all_features.output,
    output: FEATURE_TABLE_OUT_DIR / "merged_table_pivot_{antibiotic}.tsv",
    log: FEATURE_TABLES_LOGS_DIR / "pivot_merged_features_miller_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_merged_features_miller_{antibiotic}.tsv",
    conda: ENVS_DIR.format("miller"),
    container: CONTAINERS.format("miller:1.0.0")
    threads: workflow.cores,
    params:
        # The rule keeps its historical name. Miller's in-memory reshape cannot
        # hold the multi-GB merged table, so the script streams instead.
        sort_memory = lambda wc, resources: f"{max(resources.mem_mb // 2, 500)}M",
        # Sort spills about the input size; keep it on shared storage beside the output.
        sort_tmp = subpath(output[0], parent=True),
    shell:
        r"""
        # Rows are samples (field 1), columns are features (field 2), empty fill.
        bash {LONG_TO_WIDE:q} {input:q} {output:q} {params.sort_tmp:q} {threads} {params.sort_memory} \
            1 2 hash '' 2> {log:q}
        """

use rule pivot_merged_features_miller as pivot_feature_set with:
    # Not bundled: each pivot uses all cores, so a shared group job would run them in series.
    group: "pivot_feature_set_batch"
    input: rules.merge_features.output,
    output: EXPERIMENT_DIR / "feature_table" / "merged_table_pivot_{antibiotic}.tsv",
    log: EXPERIMENT_DIR / "feature_table" / "logs" / "pivot_feature_set_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_feature_set_{feature_set}_{antibiotic}.tsv",

rule feature_table:
    localrule: True
    input:
        expand(rules.pivot_merged_features_miller.output, antibiotic = ANTIBIOTICS)
    output: touch(OUT_DIR / "flags" / "feature_table.done")

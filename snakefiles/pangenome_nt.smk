from pathlib import Path

# Nucleotide pangenome experiment: SNP + GPA + nucleotide pangenome versus the
# baseline SNP + GPA + panproteome. The baseline "pangenome" rules build protein
# graphs from Prokka's translated CDS (.faa). This variant keeps the same CD-HIT
# protein clusters from the GPA step and replaces each protein with its coding
# sequence from Prokka's .ffn, so only the sequence alphabet changes. Prokka's
# .fna holds the input contigs, not genes, and CD-HIT-EST cannot cluster at the
# 0.7 identity used for proteins (its minimum is 0.75).
#
# GPA, its CD-HIT run and the baseline pangenome rules are untouched. Every rule
# here writes under pangenome_nt/ or experiments/snp_gpa_pangenome_nt/.

PANGENOME_NT_OUT_DIR = OUT_DIR / "pangenome_nt"
PANGENOME_NT_LOGS_DIR = PANGENOME_NT_OUT_DIR / "logs"
NT_SHARD_MANIFEST = PANGENOME_NT_OUT_DIR / "cluster_shards" / "{shard}.txt"
EXPERIMENT_NT_DIR = OUT_DIR / "experiments" / "snp_gpa_pangenome_nt"


# -----------------------
# Gene Sequences
# -----------------------

rule combine_ffn_files:
    localrule: True
    input:
        lambda wc: expand(
            rules.prokka_runner.output.ffn,
            sample = get_sample_names(wc)
        )
    output: PANGENOME_NT_OUT_DIR / "combined_genes.ffn",
    shell:
        # Same header rewrite as cdhit_protein_name_corrector, so CDS names match
        # the protein names in the CD-HIT clusters. The sample is the file stem.
        r"""
        awk 'BEGIN {{OFS=""}} \
            FNR == 1 {{
                sample = FILENAME;
                sub(/^.*\//, "", sample);
                sub(/\.ffn$/, "", sample);
            }} \
            !/^>/ {{print; next}} \
            /^>/ {{
                header = substr($0, 2);
                gsub(/[^a-zA-Z0-9]/, "_", header);
                sub(/_+$/, "", header);
                print ">" sample "_" header
            }}' {input} > {output}
        """


# -----------------------
# Split Clusters
# -----------------------

use rule split_cluster_fasta as split_cluster_fasta_nt with:
    input:
        cdhit_clstr = rules.cdhit_runner.output.clstr,
        # The script's input name is historical; these are nucleotide sequences.
        combined_proteins = rules.combine_ffn_files.output[0],
    output: directory(PANGENOME_NT_OUT_DIR / "cluster_sequences"),
    log: PANGENOME_NT_LOGS_DIR / "split_cluster_fasta.log",
    benchmark: BENCHMARKS_DIR / "split_cluster_fasta_nt.tsv",


use rule shard_clusters as shard_clusters_nt with:
    input: rules.split_cluster_fasta_nt.output[0],
    output: expand(NT_SHARD_MANIFEST, shard=PANGENOME_SHARDS),


use rule cluster_fasta_splits as cluster_fasta_splits_nt with:
    group: "cluster_fasta_splits_nt_batch"
    input:
        cluster_store = rules.split_cluster_fasta_nt.output[0],
        datasail_splits = rules.datasail_runner.output[0],
    output: directory(PANGENOME_NT_OUT_DIR / "cluster_fasta_splits" / "{antibiotic}" / "{split_category}"),
    log: PANGENOME_NT_LOGS_DIR / "cluster_fasta_splits" / "{antibiotic}_{split_category}.log",
    benchmark: BENCHMARKS_DIR / "cluster_fasta_splits_nt_{antibiotic}_{split_category}.tsv",


# -----------------------
# Alignment and Graphs
# -----------------------

# mafft --auto detects nucleotide input.
use rule align_clusters as align_clusters_nt with:
    input:
        manifest = NT_SHARD_MANIFEST,
        cluster_store = rules.split_cluster_fasta_nt.output[0],
    output: directory(PANGENOME_NT_OUT_DIR / "cluster_alignments" / "{shard}")
    log: PANGENOME_NT_LOGS_DIR / "align_clusters" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "align_clusters_nt_{shard}.tsv"


# PanPA builds graphs from any MSA alphabet.
use rule panpa_build_gfa as panpa_build_gfa_nt with:
    input:
        manifest = NT_SHARD_MANIFEST,
        alignments = rules.align_clusters_nt.output[0],
    output: directory(PANGENOME_NT_OUT_DIR / "panpa" / "gfa" / "{shard}")
    log: PANGENOME_NT_LOGS_DIR / "panpa_build_gfa" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "panpa_build_gfa_nt_{shard}.tsv"


use rule bubblegun_runner as bubblegun_runner_nt with:
    input:
        manifest = NT_SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa_nt.output[0],
    output: directory(PANGENOME_NT_OUT_DIR / "bubblegun" / "{shard}")
    log: PANGENOME_NT_LOGS_DIR / "bubblegun_runner" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubblegun_runner_nt_{shard}.tsv"


# -----------------------
# Bubble Features: Train
# -----------------------

# Same as bubble_features, with two options the protein run did not need.
# Nucleotide graphs have bubbles whose realized paths grow combinatorially:
# Cluster_4184 exceeded 30 GB, and debug logs reached 5.6 GB per shard
# (run 67228). Bubbles with more than BUBBLE_PATH_LIMIT paths are skipped,
# and only warnings and errors are logged.
BUBBLE_PATH_LIMIT = 5000

rule bubble_features_nt:
    input:
        manifest = NT_SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa_nt.output[0],
        bubblegun_dir = rules.bubblegun_runner_nt.output[0],
        phenotype_table = lambda wc: rules.split_phenotype_dataframe.output[0].format(
            antibiotic=wc.antibiotic, split_category="train"
        ),
    output:
        output_dir = directory(PANGENOME_NT_OUT_DIR / "bubble_features" / "train" / "{antibiotic}" / "{shard}"),
        lor_lookup_dir = directory(PANGENOME_NT_OUT_DIR / "bubble_features_{antibiotic}_lor_lookup" / "{shard}"),
    log: PANGENOME_NT_LOGS_DIR / "bubble_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubble_features_nt_{antibiotic}_{shard}.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output.output_dir:q} {output.lor_lookup_dir:q}
        export SCRIPT={BUBBLE_FEATURES_SCRIPT:q} GFA={input.gfa_dir:q} BUBBLES={input.bubblegun_dir:q} \
            PHENOTYPES={input.phenotype_table:q} ANTIBIOTIC={wildcards.antibiotic:q} \
            OUT={output.output_dir:q} LOOKUP={output.lor_lookup_dir:q} PATH_LIMIT={BUBBLE_PATH_LIMIT}
        bash {FOR_EACH_CLUSTER:q} {input.manifest:q} {threads} {log:q} '
            python "$SCRIPT" \
                --gfa-file "$GFA/$1.gfa" --bubble-gun "$BUBBLES/$1.json" \
                --phenotype-table "$PHENOTYPES" --antibiotic "$ANTIBIOTIC" --log-file "$2" \
                --output-file "$OUT/$1.tsv" --lor-lookup-file "$LOOKUP/$1.tsv" \
                --max-bubble-paths "$PATH_LIMIT" --log-level WARNING
            test -f "$OUT/$1.tsv" && test -f "$LOOKUP/$1.tsv"'
        """


use rule merge_bubble_features as merge_bubble_features_nt with:
    input: expand(rules.bubble_features_nt.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: PANGENOME_NT_OUT_DIR / "bubble_features_{antibiotic}.tsv"


# -----------------------
# Bubble Features: Test
# -----------------------

rule panpa_align_nt:
    input:
        manifest = NT_SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa_nt.output[0],
        test_sequence_folder = lambda wc: rules.cluster_fasta_splits_nt.output[0].format(
            antibiotic=wc.antibiotic, split_category="test",
        ),
    output: directory(PANGENOME_NT_OUT_DIR / "panpa" / "alignments" / "{antibiotic}" / "{shard}")
    log: PANGENOME_NT_LOGS_DIR / "panpa_align" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "panpa_align_nt_{antibiotic}_{shard}.tsv"
    params:
        # Nucleotide queries against nucleotide graphs use PanPA's DNA scoring
        # matrix (match 1, mismatch 0). PanPA's --dna flag is not used: it aligns
        # DNA queries to protein graphs by translating codons, frameshift-aware.
        sub_matrix = "dna",
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output:q}
        export GFA={input.gfa_dir:q} QUERIES={input.test_sequence_folder:q} OUT={output:q} \
            SUB_MATRIX={params.sub_matrix:q}
        bash {FOR_EACH_CLUSTER:q} {input.manifest:q} {threads} {log:q} '
            query="$QUERIES/$1" gaf="$OUT/$1.gaf"
            if [ ! -f "$query" ]; then
                echo "Missing test FASTA: $query"
                exit 1
            fi
            if [ ! -s "$query" ]; then
                echo "No test sequences; creating an empty GAF."
                : > "$gaf"
                exit 0
            fi
            PanPA --log_file "$2" align_single \
                --gfa_files "$GFA/$1.gfa" --seqs "$query" --sub_matrix "$SUB_MATRIX" \
                --cores 1 --out_gaf "$gaf"
            # PanPA can finish successfully without emitting any alignments.
            if [ ! -f "$gaf" ]; then
                : > "$gaf"
            fi'
        """


use rule gaf_lor_features as gaf_lor_features_nt with:
    input:
        manifest = NT_SHARD_MANIFEST,
        gaf_dir = rules.panpa_align_nt.output[0],
        lor_lookup_dir = rules.bubble_features_nt.output.lor_lookup_dir,
    output:
        output_dir = directory(PANGENOME_NT_OUT_DIR / "bubble_features" / "test" / "{antibiotic}" / "{shard}"),
    log: PANGENOME_NT_LOGS_DIR / "gaf_lor_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "gaf_lor_features_nt_{antibiotic}_{shard}.tsv"


use rule merge_bubble_features_test as merge_bubble_features_test_nt with:
    input: expand(rules.gaf_lor_features_nt.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: PANGENOME_NT_OUT_DIR / "bubble_features_test_{antibiotic}.tsv"


rule pangenome_nt:
    localrule: True
    input:
        expand(rules.merge_bubble_features_nt.output, antibiotic=ANTIBIOTICS),
        expand(rules.merge_bubble_features_test_nt.output, antibiotic=ANTIBIOTICS),
    output: touch(OUT_DIR / "flags" / "pangenome_nt.done")


# -----------------------
# Feature Table and ML
# -----------------------

# Same SNP and GPA tables as the baseline; only the pangenome features differ.
use rule merge_features as merge_features_nt with:
    input:
        snp = rules.binary_mutation_table.output,
        gpa = rules.binary_gpa.output,
        pangenome_train = rules.merge_bubble_features_nt.output,
        pangenome_test = rules.merge_bubble_features_test_nt.output,
    output: EXPERIMENT_NT_DIR / "feature_table" / "merged_table_{antibiotic}.tsv",


use rule pivot_merged_features_miller as pivot_features_nt with:
    group: "pivot_features_nt_batch"
    input: rules.merge_features_nt.output,
    output: EXPERIMENT_NT_DIR / "feature_table" / "merged_table_pivot_{antibiotic}.tsv",
    log: EXPERIMENT_NT_DIR / "feature_table" / "logs" / "pivot_features_nt_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_features_nt_{antibiotic}.tsv",


# Same model, grid and search as the baseline combined_ml.
EXPERIMENT_NT_ML_STEM = str(EXPERIMENT_NT_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}")

use rule combined_ml as combined_ml_nt with:
    group: "combined_ml_nt_batch"
    input:
        binary_mutation_table = rules.pivot_features_nt.output[0],
        phenotype_table = rules.phenotype_dataframe_creator.output[0],
        train = lambda wildcards: expand(rules.split_train_test.output, split_category=["train"], **wildcards)[0],
        test = lambda wildcards: expand(rules.split_train_test.output, split_category=["test"], **wildcards)[0],
    output:
        best_params = EXPERIMENT_NT_ML_STEM + "_best_params.txt",
        model_file  = EXPERIMENT_NT_ML_STEM + "_model.sav",
        result      = EXPERIMENT_NT_ML_STEM + "_Result.txt",
        fia         = EXPERIMENT_NT_ML_STEM + ".txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_nt_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.tsv"
    log: EXPERIMENT_NT_ML_STEM + ".log"


rule experiment_pangenome_nt:
    localrule: True
    input:
        expand(
            rules.combined_ml_nt.output,
            antibiotic=ANTIBIOTICS,
            random_seed=[42],
            test_size=[0.2],
            resampling_strategy=["cv"],
            model_type=["xgb"],
            feature_importance_analysis_strategy=["gini"]
        )
    output: touch(OUT_DIR / "flags" / "experiment_pangenome_nt.done")

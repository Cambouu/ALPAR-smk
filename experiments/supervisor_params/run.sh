#!/bin/bash
# One HTCondor job per antibiotic: every feature table, both splits, both parameter sets.
# Usage: run.sh <antibiotic> <outdir>
set -euo pipefail
abx=$1
outdir=$2
out=/home/joca00004/out
here=$(dirname "$(readlink -f "$0")")

# Smallest tables first so early results arrive quickly.
tables=()
for set in gpa snp snp_gpa pangenome gpa_pangenome snp_pangenome; do
    tables+=(--table "$set=$out/experiments/$set/feature_table/merged_table_pivot_$abx.tsv")
done
tables+=(--table "baseline=$out/feature_table/merged_table_pivot_$abx.tsv")
tables+=(--table "snp_gpa_panproteome_fixed=$out/experiments/snp_gpa_panproteome_fixed/feature_table/merged_table_pivot_$abx.tsv")

exec python3 "$here/fixed_params_xgb.py" \
    --antibiotic "$abx" \
    --phenotype-table "$out/data_preprocessing/phenotype_table.tsv" \
    --datasail-dir "$out/datasail" \
    --outdir "$outdir" \
    --threads 32 \
    "${tables[@]}"

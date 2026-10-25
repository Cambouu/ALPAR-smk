"""Compare the PanGraph experiment's feature sets (pangraph.smk): test metrics, feature columns,
PanGraph feature statistics, top features and cost.

Writes a Markdown report and a long TSV (set, antibiotic, metric, value). Results of earlier runs
(the baseline, experiment/feature-sets) are added as reference rows when their files exist; they
are not recomputed. Standard library only.

Feature families are recognised by name:
  pangraph  pangraph|<kind>|...                  (pangraph_features.py)
  snp       <pos>,<REF>:<ALT>[,<type>]           (binary_mutation_table.py, snippy)
  gpa       <40-hex strain checksum>_<protein>   (binary_gpa_cdhit.py, CD-HIT representative)
"""

import logging
import re
from collections import Counter, defaultdict
from contextlib import suppress
from pathlib import Path

with suppress(ImportError):
    from snakemake.script import snakemake

METRICS = {
    "Accuracy score": "accuracy",
    "Balanced Accuracy score": "balanced_accuracy",
    "Matthews correlation coefficient": "mcc",
    "F1 score binary": "f1",
    "Precision score": "precision",
    "Recall score": "recall",
    "ROC AUC Score": "roc_auc",
}
SNIPPY = re.compile(r"^(\d+),([A-Za-z]+):[^,]+(?:,\w+)?$")
GPA = re.compile(r"^[0-9a-f]{40}_")
SUMMARY = re.compile(r"strains (\d+)\s+distinct features (\d+)\s+rows (\d+)")
CARRIER_BINS = ((1, 1), (2, 2), (3, 9), (10, 99), (100, 999), (1000, None))


def parse_result(path):
    """combined_ml Result.txt -> metric -> value; confusion matrix as tn, fp, fn, tp."""
    out, lines = {}, Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        key, _, value = line.partition(": ")
        if key in METRICS:
            out[METRICS[key]] = float(value)
        elif line.startswith("Confusion matrix:"):
            cells = []
            for row in lines[i + 1:]:
                cells += [int(x) for x in re.findall(r"-?\d+", row)]
                if row.rstrip().endswith("]]"):
                    break
            if len(cells) == 4:
                out.update(zip(("tn", "fp", "fn", "tp"), cells))
    return out


def parse_best_params(path):
    for line in Path(path).read_text().splitlines():
        if line.startswith("Parameters: "):
            return line.removeprefix("Parameters: ")
    return ""


def read_importances(path, limit):
    rows = []
    for line in Path(path).read_text().splitlines():
        feature, _, value = line.rpartition("\t")
        if feature:
            rows.append((feature, float(value)))
        if len(rows) == limit:
            break
    return rows


def pivot_columns(path):
    """Feature names from the wide table's header line (the header cell is the strain column)."""
    with open(path, "rb") as handle:
        header = handle.readline().decode().rstrip("\r\n").split("\t")
    return header[1:]


def family(feature):
    if feature.startswith("pangraph|"):
        return "pangraph"
    if SNIPPY.match(feature):
        return "snp"
    if GPA.match(feature):
        return "gpa"
    return "other"


def pangraph_kind(feature):
    """blk, or <type>_ref / <type>_block (as pangraph_features.kind)."""
    _, kind, where = feature.split("|", 2)
    if kind == "blk":
        return "blk"
    return f"{kind}_block" if where.startswith("b") else f"{kind}_ref"


def reference_interval(feature):
    """1-based reference interval a feature names, or None (GPA, block-consensus variants)."""
    if match := SNIPPY.match(feature):
        pos = int(match.group(1))
        return pos, pos + len(match.group(2)) - 1
    if not feature.startswith("pangraph|"):
        return None
    parts = feature.split("|")
    kind, where = parts[1], parts[2]
    if kind == "blk":
        if len(parts) < 4:
            return None
        lo, hi = parts[3].rpartition(":")[2].split("-")
        return int(lo), int(hi)
    if where.startswith("b"):
        return None
    if kind == "snp":
        return int(where), int(where)
    if "^" in where:
        lo, hi = where.split("^")
        return int(lo), int(lo if hi == "edge" else hi)
    lo, hi = where.split("-")
    return int(lo), int(hi)


def loci_of(feature, loci):
    interval = reference_interval(feature)
    if interval is None:
        return []
    return [name for name, (start, end) in loci.items() if interval[0] <= end and interval[1] >= start]


def read_benchmark(path):
    """Last line of a Snakemake benchmark: wall seconds and peak RSS (MB), or None."""
    path = Path(path)
    if not path.is_file():
        return None
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return None
    row = dict(zip(lines[0].split("\t"), lines[-1].split("\t")))
    with suppress(KeyError, ValueError):
        rss = float(row["max_rss"]) if row.get("max_rss", "NA") != "NA" else None
        return float(row["s"]), rss
    return None


def carrier_summary(path, loci, examples=5):
    kinds, bins, rows = Counter(), Counter(), 0
    in_locus = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            feature, _, count = line.rstrip("\n").rpartition("\t")
            count = int(count)
            rows += count
            kinds[pangraph_kind(feature)] += 1
            for low, high in CARRIER_BINS:
                if count >= low and (high is None or count <= high):
                    bins[(low, high)] += 1
            for name in loci_of(feature, loci):
                in_locus[name].append((count, feature))
    top = {name: sorted(hits, key=lambda x: (-x[0], x[1]))[:examples] for name, hits in in_locus.items()}
    counts = {name: len(hits) for name, hits in in_locus.items()}
    return dict(kinds=kinds, bins=bins, rows=rows, distinct=sum(kinds.values()), locus_counts=counts, locus_top=top)


def carriers_at_least(path, minimum):
    with open(path) as handle:
        return sum(1 for line in handle if int(line.rstrip("\n").rpartition("\t")[2]) >= minimum)


def fmt(value, digits=4, sign=False):
    return "–" if value is None else f"{value:{'+' if sign else ''}.{digits}f}"


def bin_label(low, high):
    return f"≥ {low}" if high is None else (str(low) if low == high else f"{low}–{high}")


def write_report(results, best_params, fia, pivots, carriers, samples, markdown, table,
                 feature_sets, reference_set=None, min_carriers=1, features_log=None, reference_ml_dirs=None,
                 benchmarks_dir=None, loci=None, top_features=20):
    loci = {name: tuple(bounds) for name, bounds in (loci or {}).items()}
    # The reference set comes first in every table.
    set_names = list(feature_sets)
    if reference_set in feature_sets:
        set_names.remove(reference_set)
        set_names.insert(0, reference_set)

    def key(path):
        # <set>/ml/<antibiotic>/<stem>_Result.txt
        parts = Path(path).parts
        return parts[-4], parts[-2]

    runs = {key(path): {"result": path} for path in results}
    for name, paths in (("best_params", best_params), ("fia", fia)):
        for path in paths:
            runs[key(path)][name] = path
    antibiotics = sorted({antibiotic for _, antibiotic in runs})
    result_name = Path(results[0]).name
    stem = result_name.removesuffix("_Result.txt")

    metrics = {k: parse_result(v["result"]) for k, v in runs.items()}
    params = {k: parse_best_params(v["best_params"]) for k, v in runs.items()}
    references = {}
    for label, directory in (reference_ml_dirs or {}).items():
        for antibiotic in antibiotics:
            path = Path(directory) / antibiotic / result_name
            if path.is_file():
                references[label, antibiotic] = parse_result(path)
            else:
                logging.info("reference result not found: %s", path)

    columns = {Path(path).parts[-3]: pivot_columns(path) for path in pivots}
    column_families = {name: Counter(map(family, cols)) for name, cols in columns.items()}
    column_kinds = {name: Counter(pangraph_kind(c) for c in cols if c.startswith("pangraph|"))
                    for name, cols in columns.items()}

    n_samples = sum(1 for line in Path(samples).read_text().splitlines() if line.strip())
    stats = carrier_summary(carriers, loci)
    kept = stats["distinct"] if min_carriers <= 1 else carriers_at_least(carriers, min_carriers)
    graph_strains = None
    if features_log and Path(features_log).is_file():
        if match := SUMMARY.search(Path(features_log).read_text()):
            graph_strains = int(match.group(1))

    tsv = []
    md = []
    title = " vs ".join(" + ".join(f.upper() if f != "pangraph" else "PanGraph" for f in feature_sets[s])
                        for s in set_names)
    md.append(f"# PanGraph experiment: {title}\n")
    md.append(
        "Test split of each antibiotic's DataSAIL split; XGBoost, small grid (36 combinations), seed 42, "
        "low-memory search forced for every set (hyperparameters are selected on the test split, as in "
        f"the baseline). PanGraph features: minimum carriers {min_carriers} (label-free). Model files: `{stem}_*`.\n"
    )

    md.append("## Test-split results\n")
    md.append("| Feature set | Antibiotic | Balanced accuracy | MCC | Accuracy | F1 | ROC AUC | TN | FP | FN | TP | Best parameters |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")

    def metric_row(label, antibiotic, m, best=""):
        cells = [fmt(m.get(x)) for x in ("balanced_accuracy", "mcc", "accuracy", "f1", "roc_auc")]
        cells += [str(m.get(x, "–")) for x in ("tn", "fp", "fn", "tp")]
        md.append(f"| {label} | {antibiotic} | " + " | ".join(cells) + f" | {best} |")

    for antibiotic in antibiotics:
        for name in set_names:
            metric_row(f"**{name}**", antibiotic, metrics[name, antibiotic], params[name, antibiotic])
            for metric, value in metrics[name, antibiotic].items():
                tsv.append((name, antibiotic, metric, value))
        for (label, ab), m in references.items():
            if ab == antibiotic:
                metric_row(f"{label} (earlier run)", antibiotic, m)
                for metric, value in m.items():
                    tsv.append((label, antibiotic, metric, value))
    md.append("")

    if len(set_names) > 1:
        base = set_names[0]
        md.append(f"## Differences against {base}\n")
        md.append("Positive Δ means the set scores higher. With about 180 (ampicillin) and 314 (ciprofloxacin) "
                  "test strains and one seed, a difference of one or two strains is not evidence of a real difference.\n")
        md.append("| Antibiotic | Feature set | Δ balanced accuracy | Δ MCC | Δ FP | Δ FN |")
        md.append("|---|---|---|---|---|---|")
        for antibiotic in antibiotics:
            reference = metrics[base, antibiotic]
            for name in set_names[1:]:
                m = metrics[name, antibiotic]
                delta = {x: (m[x] - reference[x]) if x in m and x in reference else None
                         for x in ("balanced_accuracy", "mcc", "fp", "fn")}
                md.append(f"| {antibiotic} | {name} | {fmt(delta['balanced_accuracy'], sign=True)} | "
                          f"{fmt(delta['mcc'], sign=True)} | {fmt(delta['fp'], 0, True)} | {fmt(delta['fn'], 0, True)} |")
        md.append("")

    md.append("## Feature columns\n")
    families = ["snp", "gpa", "pangraph", "other"]
    md.append("| Feature set | Columns | " + " | ".join(families) + " |")
    md.append("|---|---|" + "---|" * len(families))
    for name in set_names:
        counts = column_families.get(name, Counter())
        md.append(f"| {name} | {sum(counts.values()):,} | " + " | ".join(f"{counts[f]:,}" for f in families) + " |")
        for f in families:
            tsv.append((name, "", f"columns_{f}", counts[f]))
    for name, kinds in column_kinds.items():
        if kinds:
            md.append(f"\nPanGraph columns in {name}: " + ", ".join(f"{k} {v:,}" for k, v in sorted(kinds.items())) + ".")
    md.append("")

    md.append("## PanGraph features\n")
    strains_line = f"{n_samples:,} strains in the graph input"
    if graph_strains is not None:
        strains_line += f"; {graph_strains:,} with features"
        if graph_strains != n_samples:
            strains_line += " (**mismatch**)"
    md.append(f"- {strains_line}.")
    md.append(f"- {stats['distinct']:,} distinct features, {stats['rows']:,} strain-feature rows before the filter; "
              f"{kept:,} features carried by at least {min_carriers} strain(s) are in the tables.")
    md.append("- By kind: " + ", ".join(f"{k} {v:,}" for k, v in sorted(stats["kinds"].items())) + ".")
    md.append("- Carriers per feature: " + ", ".join(
        f"{bin_label(*b)}: {stats['bins'][b]:,}" for b in CARRIER_BINS) + ".")
    for name in loci:
        top = stats["locus_top"].get(name, [])
        md.append(f"- {name} ({loci[name][0]:,}–{loci[name][1]:,}): {stats['locus_counts'].get(name, 0):,} features; most carried: "
                  + (", ".join(f"`{f}` ({c:,})" for c, f in top) or "none") + ".")
    tsv += [("pangraph", "", "graph_input_strains", n_samples), ("pangraph", "", "distinct_features", stats["distinct"]),
            ("pangraph", "", "feature_rows", stats["rows"]), ("pangraph", "", f"features_min_carriers_{min_carriers}", kept)]
    if graph_strains is not None:
        tsv.append(("pangraph", "", "strains_with_features", graph_strains))
    md.append("")

    md.append(f"## Top {top_features} features (XGBoost weight)\n")
    for antibiotic in antibiotics:
        for name in set_names:
            rows = read_importances(runs[name, antibiotic]["fia"], top_features)
            hits = sum(1 for feature, _ in rows if loci_of(feature, loci))
            md.append(f"**{name} / {antibiotic}** ({hits} of {len(rows)} in " + ", ".join(loci) + ")\n")
            md.append("| Rank | Feature | Weight | Family | Locus |")
            md.append("|---|---|---|---|---|")
            for rank, (feature, value) in enumerate(rows, 1):
                md.append(f"| {rank} | `{feature}` | {value:g} | {family(feature)} | {', '.join(loci_of(feature, loci))} |")
            md.append("")
            tsv.append((name, antibiotic, f"top{top_features}_in_loci", hits))

    md.append("## Cost (Snakemake benchmarks)\n")
    md.append("| Step | Wall time (h) | Peak RSS (GiB) |")
    md.append("|---|---|---|")
    steps = [(rule, f"{rule}.tsv") for rule in ("pangraph_input", "pangraph_build", "pangraph_features")]
    steps += [(f"pivot {name}", f"pivot_pangraph_set_{name}.tsv") for name in set_names]
    steps += [(f"combined_ml {name} / {antibiotic}", f"combined_ml_pangraph_{name}_{antibiotic}_{stem}.tsv")
              for antibiotic in antibiotics for name in set_names]
    for label, filename in steps:
        bench = read_benchmark(Path(benchmarks_dir) / filename) if benchmarks_dir else None
        if bench is None:
            md.append(f"| {label} | – | – |")
            continue
        seconds, rss = bench
        md.append(f"| {label} | {seconds / 3600:.2f} | {fmt(None if rss is None else rss / 1024, 1)} |")
        tsv.append((label, "", "wall_s", seconds))
        if rss is not None:
            tsv.append((label, "", "max_rss_mb", rss))
    md.append("\nHTCondor jobs record peak memory more reliably (`condor_history -af MemoryUsage`); "
              "benchmarks of jobs run in containers have reported 0.\n")

    Path(markdown).write_text("\n".join(md) + "\n", encoding="utf-8")
    with open(table, "w", encoding="utf-8") as handle:
        handle.write("set\tantibiotic\tmetric\tvalue\n")
        for row in tsv:
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)
    write_report(
        results=list(snakemake.input.results),
        best_params=list(snakemake.input.best_params),
        fia=list(snakemake.input.fia),
        pivots=list(snakemake.input.pivots),
        carriers=snakemake.input.carriers,
        samples=snakemake.input.samples,
        markdown=snakemake.output.markdown,
        table=snakemake.output.table,
        feature_sets=snakemake.params.feature_sets,
        reference_set=snakemake.params.reference_set,
        min_carriers=int(snakemake.params.min_carriers),
        features_log=snakemake.params.features_log,
        reference_ml_dirs=snakemake.params.reference_ml_dirs,
        benchmarks_dir=snakemake.params.benchmarks_dir,
        loci=snakemake.params.loci,
        top_features=int(snakemake.params.top_features),
    )

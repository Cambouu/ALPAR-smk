"""Report of the Minigraph-Cactus experiment (minigraph_cactus.smk): test metrics against the
reference runs of experiment/pangraph, feature columns, MC feature statistics, sequence kept in the
graph, top features and cost.

Writes a Markdown report and a long TSV (set, antibiotic, metric, value). Standard library only.

Feature families are recognised by name:
  mc        mc|<snp|ind>|<pos>|<REF>><ALT>       (mc_vcf_to_features.py)
  pangraph  pangraph|<kind>|...                  (pangraph_features.py)
  snp       <pos>,<REF>:<ALT>[,<type>]           (binary_mutation_table.py, snippy)
  gpa       <40-hex strain checksum>_<protein>   (binary_gpa_cdhit.py, CD-HIT representative)
"""

import logging
import re
import statistics
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
FAMILY_LABELS = {"gpa": "GPA", "mc": "MC", "snp": "SNP", "pangraph": "PanGraph"}


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
            with suppress(ValueError):
                rows.append((feature, float(value)))
        if len(rows) == limit:
            break
    return rows


def pivot_columns(path):
    with open(path, "rb") as handle:
        return handle.readline().decode().rstrip("\r\n").split("\t")[1:]


def family(feature):
    if feature.startswith("mc|"):
        return "mc"
    if feature.startswith("pangraph|"):
        return "pangraph"
    if GPA.match(feature):
        return "gpa"
    if SNIPPY.match(feature):
        return "snp"
    return "other"


def position(feature):
    """Reference position of an MC or snippy feature, else None."""
    if feature.startswith("mc|"):
        with suppress(IndexError, ValueError):
            return int(feature.split("|")[2])
    m = SNIPPY.match(feature)
    return int(m.group(1)) if m else None


def loci_of(feature, loci):
    pos = position(feature)
    return ",".join(name for name, (lo, hi) in loci.items() if pos is not None and lo <= pos <= hi)


def read_benchmark(path):
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


def fmt(value, digits=4, sign=False):
    if value is None:
        return "–"
    return f"{value:+.{digits}f}" if sign else f"{value:.{digits}f}"


def feature_stats(path):
    distinct, strains, rows, kinds = set(), set(), 0, Counter()
    with open(path) as handle:
        for line in handle:
            strain, feature, _ = line.rstrip("\n").split("\t")
            rows += 1
            strains.add(strain)
            if feature not in distinct:
                distinct.add(feature)
                kinds[feature.split("|")[1]] += 1
    return {"distinct": len(distinct), "rows": rows, "strains": len(strains), "kinds": kinds}


def sequence_kept(seqfile, pathlen, ref_name):
    """Fraction of each genome's bases on its paths in the graph (vg paths -E, PanSN names)."""
    kept = Counter()
    for line in open(pathlen):
        name, _, length = line.rstrip("\n").rpartition("\t")
        kept[name.split("#")[0]] += int(length)
    fractions = {}
    for line in open(seqfile):
        sample, _, fasta = line.rstrip("\n").partition("\t")
        if sample == ref_name or not fasta:
            continue
        total = 0
        with open(fasta) as handle:
            for row in handle:
                if not row.startswith(">"):
                    total += len(row.strip())
        if total:
            fractions[sample] = min(kept.get(sample, 0) / total, 1.0)
    return fractions


def write_report(results, best_params, fia, pivots, features, samples, pathlen, seqfile, markdown, table,
                 feature_sets, vcf="wave", features_log=None, reference_ml_dirs=None, benchmarks_dir=None,
                 loci=None, top_features=20, ref_name="MG1655"):
    loci = {name: tuple(bounds) for name, bounds in (loci or {}).items()}
    set_names = list(feature_sets)

    def key(path):
        parts = Path(path).parts  # <set>/ml/<antibiotic>/<stem>_Result.txt
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

    md, tsv = [], []
    title = ", ".join(" + ".join(FAMILY_LABELS.get(f, f) for f in feature_sets[s]) for s in set_names)
    md.append(f"# Minigraph-Cactus experiment: {title}\n")
    md.append("Test split of each antibiotic's DataSAIL split; XGBoost, small grid (36 combinations), seed 42, "
              "low-memory search (hyperparameters selected on the test split, as in experiment/pangraph). "
              f"MC features from the `{vcf}` deconstruct VCF. Model files: `{stem}_*`.\n")

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
            tsv += [(name, antibiotic, metric, value) for metric, value in metrics[name, antibiotic].items()]
        for (label, ab), m in references.items():
            if ab == antibiotic:
                metric_row(f"{label} (earlier run)", antibiotic, m)
                tsv += [(label, antibiotic, metric, value) for metric, value in m.items()]
    md.append("")

    first_ref = next(iter(reference_ml_dirs or {}), None)
    if first_ref and any(label == first_ref for label, _ in references):
        md.append(f"## Differences against {first_ref}\n")
        md.append("Positive Δ means the set scores higher. With about 180 (ampicillin) and 314 (ciprofloxacin) "
                  "test strains and one seed, a difference of one or two strains is not evidence of a real difference.\n")
        md.append("| Antibiotic | Feature set | Δ balanced accuracy | Δ MCC | Δ FP | Δ FN |")
        md.append("|---|---|---|---|---|---|")
        for antibiotic in antibiotics:
            reference = references.get((first_ref, antibiotic))
            if not reference:
                continue
            for name in set_names:
                m = metrics[name, antibiotic]
                d = {x: (m[x] - reference[x]) if x in m and x in reference else None
                     for x in ("balanced_accuracy", "mcc", "fp", "fn")}
                md.append(f"| {antibiotic} | {name} | {fmt(d['balanced_accuracy'], sign=True)} | "
                          f"{fmt(d['mcc'], sign=True)} | {fmt(d['fp'], 0, True)} | {fmt(d['fn'], 0, True)} |")
        md.append("")

    md.append("## Feature columns\n")
    families = ["mc", "gpa", "snp", "pangraph", "other"]
    md.append("| Feature set | Columns | " + " | ".join(FAMILY_LABELS.get(f, f) for f in families) + " |")
    md.append("|---|---|" + "---|" * len(families))
    for path in pivots:
        name = Path(path).parts[-3]
        counts = Counter(family(c) for c in pivot_columns(path))
        md.append(f"| {name} | {sum(counts.values()):,} | " + " | ".join(f"{counts[f]:,}" for f in families) + " |")
        tsv += [(name, "", f"columns_{f}", counts[f]) for f in families]
    md.append("")

    stats = feature_stats(features)
    n_samples = sum(1 for line in open(samples) if line.strip())
    md.append("## MC features\n")
    line = f"- {n_samples:,} strains in the graph input; {stats['strains']:,} with features"
    md.append(line + (" (**mismatch**)." if stats["strains"] != n_samples else "."))
    md.append(f"- {stats['distinct']:,} distinct features ({', '.join(f'{k} {v:,}' for k, v in sorted(stats['kinds'].items()))}), "
              f"{stats['rows']:,} strain-feature rows.")
    if features_log and Path(features_log).is_file():
        md.append("- Extraction log: `" + Path(features_log).read_text().strip().splitlines()[-1] + "`.")
    tsv += [("mc", "", "graph_input_strains", n_samples), ("mc", "", "strains_with_features", stats["strains"]),
            ("mc", "", "distinct_features", stats["distinct"]), ("mc", "", "feature_rows", stats["rows"])]
    tsv += [("mc", "", f"distinct_{k}", v) for k, v in stats["kinds"].items()]
    fractions = sequence_kept(seqfile, pathlen, ref_name)
    if fractions:
        values = sorted(fractions.values())
        md.append(f"- Sequence kept in the graph per genome: median {statistics.median(values):.1%}, "
                  f"range {values[0]:.1%}–{values[-1]:.1%}; {sum(v < 0.8 for v in values)} genomes below 80%.")
        tsv += [("mc", "", "kept_median", statistics.median(values)), ("mc", "", "kept_min", values[0]),
                ("mc", "", "kept_max", values[-1])]
    md.append("")

    md.append(f"## Top {top_features} features\n")
    for antibiotic in antibiotics:
        for name in set_names:
            path = runs[name, antibiotic].get("fia")
            if not path:
                continue
            md.append(f"**{name} / {antibiotic}**\n")
            md.append("| Rank | Feature | Family | Locus | Importance |")
            md.append("|---|---|---|---|---|")
            for rank, (feature, value) in enumerate(read_importances(path, top_features), 1):
                md.append(f"| {rank} | `{feature}` | {family(feature)} | {loci_of(feature, loci) or '–'} | {value:.4f} |")
            md.append("")

    md.append("## Cost (Snakemake benchmarks)\n")
    md.append("| Step | Wall time (h) | Peak RSS (GB) |")
    md.append("|---|---|---|")
    benchmarks = Path(benchmarks_dir) if benchmarks_dir else None

    def row(label, file_name):
        bench = read_benchmark(benchmarks / file_name) if benchmarks else None
        if bench is None:
            md.append(f"| {label} | – | – |")
            return
        wall, rss = bench
        md.append(f"| {label} | {wall / 3600:.2f} | {'–' if rss is None else f'{rss / 1000:.1f}'} |")
        tsv.append((label, "", "wall_s", wall))
        if rss is not None:
            tsv.append((label, "", "max_rss_mb", rss))

    row("mc_input", "mc_input.tsv")
    row("mc_build", "mc_build.tsv")
    row("mc_features", "mc_features.tsv")
    for name in set_names:
        row(f"pivot {name}", f"pivot_mc_set_{name}.tsv")
        for antibiotic in antibiotics:
            row(f"combined_ml {name} / {antibiotic}", f"combined_ml_mc_{name}_{antibiotic}_{stem}.tsv")
    md.append("")

    Path(markdown).write_text("\n".join(md) + "\n")
    with open(table, "w") as out:
        out.write("set\tantibiotic\tmetric\tvalue\n")
        for record in tsv:
            out.write("\t".join(str(x) for x in record) + "\n")


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)
    write_report(
        results=list(snakemake.input.results),
        best_params=list(snakemake.input.best_params),
        fia=list(snakemake.input.fia),
        pivots=list(snakemake.input.pivots),
        features=snakemake.input.features,
        samples=snakemake.input.samples,
        pathlen=snakemake.input.pathlen,
        seqfile=snakemake.input.seqfile,
        markdown=snakemake.output.markdown,
        table=snakemake.output.table,
        feature_sets=snakemake.params.feature_sets,
        vcf=snakemake.params.vcf,
        features_log=snakemake.params.features_log,
        reference_ml_dirs=snakemake.params.reference_ml_dirs,
        benchmarks_dir=snakemake.params.benchmarks_dir,
        loci=snakemake.params.loci,
        top_features=int(snakemake.params.top_features),
    )

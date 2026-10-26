"""Tests for pangraph_report.py on small hand-made ML outputs (standard library only)."""

from pathlib import Path
from tempfile import TemporaryDirectory
import importlib.util
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "snakefiles" / "scripts" / "pangraph_report.py"
SPEC = importlib.util.spec_from_file_location("pangraph_report", SCRIPT)
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)

STEM = "seed_42_testsize_0.2_resampling_cv_xgb_FIA_gini"
GPA = "0123456789abcdef0123456789abcdef01234567_ABCDEFGH_00001_DNA_gyrase_subunit_A"
LOCI = {"gyrA": [2336793, 2339420], "parC": [3163715, 3165973]}


def result_text(ba, mcc, tn, fp, fn, tp):
    return (
        f"Accuracy score: 0.9\nBalanced Accuracy score: {ba}\nMatthews correlation coefficient: {mcc}\n"
        "F1 score binary: 0.8\nPrecision score: 0.8\nRecall score: 0.8\n"
        f"Confusion matrix:\n[[{tn:3d} {fp:3d}]\n [{fn:3d} {tp:3d}]]\nROC AUC Score: {ba}\nBrier score loss: 0.1\n"
    )


class PangraphReportTest(unittest.TestCase):
    def build(self, root):
        sets = {"gpa_pangraph": ["gpa", "pangraph"], "snp_gpa": ["snp", "gpa"]}
        tables = {
            "snp_gpa": (result_text(0.9, 0.8, 210, 9, 6, 89), ["2339166,G:A,snp", GPA, "100,A:AT,ins"]),
            "gpa_pangraph": (result_text(0.95, 0.85, 212, 7, 7, 88),
                             [GPA, "pangraph|snp|2339173|G>A", "pangraph|blk|b5|NC_000913.3:3164000-3166000",
                              "pangraph|ins|b7:4|AC"]),
        }
        files = {"results": [], "best_params": [], "fia": [], "pivots": []}
        for name, (result, columns) in tables.items():
            ml = root / "experiments" / "pangraph" / name / "ml" / "drug"
            ml.mkdir(parents=True)
            (ml / f"{STEM}_Result.txt").write_text(result)
            (ml / f"{STEM}_best_params.txt").write_text(
                "Best MCC result for drug: 0.8\nParameters: max_depth=3, eta=0.1\n")
            (ml / f"{STEM}.txt").write_text("".join(f"{c}\t{len(columns) - i}.0\n" for i, c in enumerate(columns)))
            pivot = root / "experiments" / "pangraph" / name / "feature_table" / "merged_table_pivot.tsv"
            pivot.parent.mkdir(parents=True)
            pivot.write_text("hash\t" + "\t".join(columns) + "\nstrain1\t" + "\t".join("1" for _ in columns) + "\n")
            files["results"].append(str(ml / f"{STEM}_Result.txt"))
            files["best_params"].append(str(ml / f"{STEM}_best_params.txt"))
            files["fia"].append(str(ml / f"{STEM}.txt"))
            files["pivots"].append(str(pivot))
        carriers = root / "feature_carriers.tsv"
        carriers.write_text(
            "pangraph|blk|b5|NC_000913.3:3164000-3166000\t3\npangraph|ins|b7:4|AC\t1\n"
            "pangraph|snp|2339173|G>A\t40\npangraph|snp|5|A>C\t2\n")
        samples = root / "samples.txt"
        samples.write_text("a\nb\nc\n")
        log = root / "features.log"
        log.write_text("strains 2  distinct features 4  rows 46\n")
        baseline = root / "ml" / "drug"
        baseline.mkdir(parents=True)
        (baseline / f"{STEM}_Result.txt").write_text(result_text(0.97, 0.95, 217, 2, 5, 90))
        benchmarks = root / "benchmarks"
        benchmarks.mkdir()
        for name, seconds, rss in (("subgraph_0", 3600, 2048), ("subgraph_1", 7200, 1024),
                                   ("merge_0-1", 1800, 4096), ("merge_0-3", 900, 512), ("merge_2-3", 600, 256)):
            (benchmarks / f"pangraph_{name}.tsv").write_text(f"s\th:m:s\tmax_rss\n{seconds}.0\t0\t{rss}.0\n")
        out = {"markdown": root / "report.md", "table": root / "report.tsv"}
        report.write_report(
            **files, carriers=str(carriers), samples=str(samples), **out, feature_sets=sets,
            reference_set="snp_gpa", min_carriers=2,
            features_log=str(log), benchmarks_dir=str(benchmarks), loci=LOCI, top_features=3,
            reference_ml_dirs={"baseline": str(root / "ml"), "missing": str(root / "nowhere")},
        )
        rows = [line.split("\t") for line in out["table"].read_text(encoding="utf-8").splitlines()[1:]]
        return out["markdown"].read_text(encoding="utf-8"), {(s, a, m): v for s, a, m, v in rows}

    def test_metrics_differences_and_references(self):
        with TemporaryDirectory() as temp:
            markdown, table = self.build(Path(temp))
        self.assertEqual(table["gpa_pangraph", "drug", "mcc"], "0.85")
        self.assertEqual(table["snp_gpa", "drug", "fn"], "6")
        self.assertEqual(table["baseline", "drug", "tp"], "90")
        self.assertNotIn(("missing", "drug", "mcc"), table)
        self.assertIn("# PanGraph experiment: SNP + GPA vs GPA + PanGraph", markdown)
        self.assertIn("| drug | gpa_pangraph | +0.0500 | +0.0500 | -2 | +1 |", markdown)
        self.assertIn("| baseline (earlier run) | drug |", markdown)
        self.assertIn("max_depth=3, eta=0.1", markdown)

    def test_feature_families_carriers_and_loci(self):
        with TemporaryDirectory() as temp:
            markdown, table = self.build(Path(temp))
        self.assertEqual(table["snp_gpa", "", "columns_snp"], "2")
        self.assertEqual(table["snp_gpa", "", "columns_gpa"], "1")
        self.assertEqual(table["gpa_pangraph", "", "columns_pangraph"], "3")
        self.assertEqual(table["pangraph", "", "feature_rows"], "46")
        self.assertEqual(table["pangraph", "", "features_min_carriers_2"], "3")
        # Top 3: snippy's gyrA SNP; PanGraph's gyrA SNP and a parC block (the block-column insertion has no coordinate).
        self.assertEqual(table["snp_gpa", "drug", "top3_in_loci"], "1")
        self.assertEqual(table["gpa_pangraph", "drug", "top3_in_loci"], "2")
        self.assertIn("2 with features (**mismatch**)", markdown)
        self.assertIn("gyrA (2,336,793–2,339,420): 1 features; most carried: `pangraph|snp|2339173|G>A` (40)", markdown)
        self.assertIn("| pangraph_subgraph ×2 (longest; total 3.00 h) | 2.00 | 2.0 |", markdown)
        # Merges in tree order, then extraction (no benchmark here).
        self.assertIn("| pangraph_merge 0-1 | 0.50 | 4.0 |\n| pangraph_merge 2-3 | 0.17 | 0.2 |\n"
                      "| pangraph_merge 0-3 | 0.25 | 0.5 |\n| pangraph_features | – | – |", markdown)
        self.assertEqual(table["pangraph_subgraph", "", "total_wall_s"], "10800.0")

    def test_reference_intervals(self):
        self.assertEqual(report.reference_interval("2339166,G:A,snp"), (2339166, 2339166))
        self.assertEqual(report.reference_interval("100,ACG:A,del"), (100, 102))
        self.assertEqual(report.reference_interval("pangraph|del|6-7"), (6, 7))
        self.assertEqual(report.reference_interval("pangraph|ins|9^10|GG"), (9, 10))
        self.assertEqual(report.reference_interval("pangraph|ins|26^edge|A"), (26, 26))
        self.assertEqual(report.reference_interval("pangraph|blk|b1|chr:11-20"), (11, 20))
        self.assertIsNone(report.reference_interval("pangraph|snp|b400:3|G>A"))
        self.assertIsNone(report.reference_interval("pangraph|blk|b400"))
        self.assertIsNone(report.reference_interval(GPA))


if __name__ == "__main__":
    unittest.main()

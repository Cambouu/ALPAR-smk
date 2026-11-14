"""Unit tests for the Minigraph-Cactus experiment scripts: feature naming, genotype parsing and the
report on hand-made inputs. The bcftools pipeline itself runs only if bcftools is on PATH."""

from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import shutil
import subprocess
import sys
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "snakefiles" / "scripts"))

import mc_report  # noqa: E402
import mc_vcf_to_features as features  # noqa: E402

A = "a" * 40
B = "b" * 40


class FeatureNameTest(unittest.TestCase):
    def test_kinds_and_long_alleles(self):
        self.assertEqual(features.feature_name("2339172", "C", "T"), "mc|snp|2339172|C>T")
        self.assertEqual(features.feature_name("10", "CA", "C"), "mc|ind|10|CA>C")
        long = "A" * 51
        name = features.feature_name("5", "A", long)
        self.assertRegex(name, r"^mc\|ind\|5\|A>len51_[0-9a-f]{8}$")
        for bad in "[]< \t":
            self.assertNotIn(bad, name)

    def test_carriers_take_any_alt_allele_and_skip_missing(self):
        gts = [f"{A}=1", f"{B}=.", "c=0", "d=0|1", "e=1/1", "f=2"]
        self.assertEqual(list(features.carriers(gts)), [A, "d", "e"])


VCF = """##fileformat=VCFv4.2
##contig=<ID=NC_000913,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a}\t{b}\t_MINIGRAPH_
NC_000913\t3\t.\tGT\tAC,G\t60\t.\t.\tGT\t1\t2\t1
NC_000913\t10\t.\tC\tA\t60\t.\t.\tGT\t0\t.\t1
"""


@unittest.skipUnless(shutil.which("bcftools"), "bcftools not on PATH")
class VcfToFeaturesTest(unittest.TestCase):
    def test_split_atomised_rows_per_strain(self):
        with TemporaryDirectory() as temp:
            temp = Path(temp)
            (temp / "ref.fa").write_text(">NC_000913\nACGTACGTACGTACGTACGT\n")
            vcf = temp / "in.vcf"
            vcf.write_text(VCF.format(a=A, b=B))
            subprocess.run(["bgzip", "-f", str(vcf)], check=True)
            subprocess.run(["tabix", "-p", "vcf", str(vcf) + ".gz"], check=True)
            out = temp / "features.tsv"
            features.main([str(vcf) + ".gz", str(temp / "ref.fa"), str(out)])
            rows = [line.split("\t") for line in out.read_text().splitlines()]
            # GT>AC atomises to G>A plus T>C; strain B carries the deletion GT>G; _MINIGRAPH_ is dropped.
            self.assertEqual(sorted((s, f) for s, f, _ in rows), sorted([
                (A, "mc|snp|3|G>A"), (A, "mc|snp|4|T>C"), (B, "mc|ind|3|GT>G"),
            ]))


class ReportTest(unittest.TestCase):
    def test_report_on_small_inputs(self):
        with TemporaryDirectory() as temp:
            temp = Path(temp)
            stem = "seed_42_testsize_0.2_resampling_cv_xgb_FIA_gini"
            ml = temp / "gpa_mc" / "ml" / "drug"
            ml.mkdir(parents=True)
            result = (
                "Accuracy score: 0.9\nBalanced Accuracy score: 0.85\nMatthews correlation coefficient: 0.7\n"
                "Confusion matrix:\n[[8 1]\n [1 10]]\n"
            )
            (ml / f"{stem}_Result.txt").write_text(result)
            (ml / f"{stem}_best_params.txt").write_text("Parameters: {'max_depth': 3}\n")
            (ml / f"{stem}.txt").write_text(f"mc|snp|2339172|C>T\t0.5\n{A}_protein\t0.2\n")
            ref = temp / "ref" / "drug"
            ref.mkdir(parents=True)
            (ref / f"{stem}_Result.txt").write_text(result.replace("0.85", "0.80"))
            pivot = temp / "gpa_mc" / "feature_table" / "merged_table_pivot.tsv"
            pivot.parent.mkdir(parents=True)
            pivot.write_text(f"strain\tmc|snp|2339172|C>T\tmc|ind|5|CA>C\t{A}_protein\n")
            feats = temp / "features.tsv"
            feats.write_text(f"{A}\tmc|snp|2339172|C>T\t1\n{B}\tmc|ind|5|CA>C\t1\n")
            samples = temp / "samples.txt"
            samples.write_text(f"/x/{A}\n/x/{B}\n")
            genome = temp / "a.fa"
            genome.write_text(">c\nACGTACGTAC\n")
            seqfile = temp / "seqfile.txt"
            seqfile.write_text(f"MG1655\t/none\n{A}\t{genome}\n")
            pathlen = temp / "pathlen.tsv"
            pathlen.write_text(f"{A}#0#c#0\t8\nMG1655#0#NC_000913\t20\n")
            md, tsv = temp / "report.md", temp / "report.tsv"
            mc_report.write_report(
                results=[str(ml / f"{stem}_Result.txt")], best_params=[str(ml / f"{stem}_best_params.txt")],
                fia=[str(ml / f"{stem}.txt")], pivots=[str(pivot)], features=str(feats), samples=str(samples),
                pathlen=str(pathlen), seqfile=str(seqfile), markdown=str(md), table=str(tsv),
                feature_sets={"gpa_mc": ["gpa", "mc"]}, reference_ml_dirs={"SNP + GPA": str(temp / "ref")},
                benchmarks_dir=str(temp), loci={"gyrA": [2336793, 2339420]},
            )
            text = md.read_text()
            self.assertIn("| **gpa_mc** | drug | 0.8500 | 0.7000 |", text)
            self.assertIn("| drug | gpa_mc | +0.0500 |", text)
            self.assertIn("| gpa_mc | 3 | 2 | 1 |", text)
            self.assertIn("median 80.0%", text)
            self.assertIn("| 1 | `mc|snp|2339172|C>T` | mc | gyrA |", text)
            self.assertIn("2 strains in the graph input; 2 with features.", text)
            self.assertIn("gpa_mc\tdrug\tbalanced_accuracy\t0.85", tsv.read_text())


if __name__ == "__main__":
    unittest.main()

"""Checks for the Minigraph-Cactus experiment rules (minigraph_cactus.smk); no tools or jobs are run."""

from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import os
import subprocess
import sys
import unittest

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ConfigSettings, ResourceSettings
from snakemake.utils import format as format_command


REPO = Path(__file__).resolve().parents[1]
SAMPLE = "0123456789abcdef0123456789abcdef01234567"
OTHER = "fedcba9876543210fedcba9876543210fedcba98"


class MinigraphCactusWorkflowTest(unittest.TestCase):
    @contextmanager
    def workflow(self, **overrides):
        with TemporaryDirectory(prefix="alpar-mc-") as temp:
            root = Path(temp)
            (root / "input" / "drug").mkdir(parents=True)
            with SnakemakeApi() as api:
                workflow_api = api.workflow(
                    ResourceSettings(cores=2),
                    config_settings=ConfigSettings(config={
                        "input_dir": str(root / "input"),
                        "output_dir": str(root / "output"),
                        "env_dir": None,
                        **overrides,
                    }),
                    snakefile=REPO / "Snakefile",
                    workdir=root,
                )
                yield root, {rule.name: rule for rule in workflow_api._workflow.rules}

    def test_default_set_merges_gpa_with_the_mc_features(self):
        with self.workflow() as (_, rules):
            merge = rules["merge_mc_set"].input[0]
            self.assertEqual(merge(SimpleNamespace(mc_set="gpa_mc")), [
                rules["binary_gpa"].output[0], rules["mc_features"].output[0],
            ])
            self.assertEqual(rules["pivot_mc_set"].shellcmd, rules["pivot_merged_features_miller"].shellcmd)
            # One set, one antibiotic in the test input.
            self.assertEqual(len(rules["mc_report"].input.results), 1)

    def test_configured_vcf_and_sets_reach_the_rules(self):
        with self.workflow(mc_vcf="clip", mc_feature_sets={"snp_mc": ["snp", "mc"]}) as (_, rules):
            self.assertEqual(Path(rules["mc_features"].input.vcf).name, "mc.vcf.gz")
            self.assertIn("snp_mc", rules["mc_report"].params.feature_sets)
            self.assertEqual(rules["mc_report"].params.vcf, "clip")
        with self.workflow() as (_, rules):
            self.assertEqual(Path(rules["mc_features"].input.vcf).name, "mc.full.vcf.gz")

    def test_build_is_one_large_job_on_the_pinned_image(self):
        with self.workflow() as (_, rules):
            build = rules["mc_build"]
            self.assertIsNone(build.conda_env)
            self.assertTrue(build.container_img.startswith(
                "docker://docker.io/cambouu/alpar-smk-minigraph-cactus@sha256:"))
            self.assertEqual(build.resources["_cores"], 64)
            self.assertEqual(build.resources["mem_mb"], 512000)
            self.assertEqual(build.resources["htcondor_request_mem_mb"], 512000 * 1000 * 1000 // 1024 ** 2)
            self.assertEqual(build.params.extra, "--noSplit")
            for text in ("--vcf clip full", "--reference {params.ref_name:q}", 'minGAFBlockLength="0"',
                         "--maxCores {threads}", "/usr/local/cactus/bin"):
                self.assertIn(text, build.shellcmd)
            self.assertIn("/usr/local/cactus/bin", rules["mc_features"].shellcmd)

    def test_ml_matches_the_pangraph_experiment(self):
        with self.workflow() as (_, rules):
            params = rules["combined_ml_mc"].params
            self.assertTrue(params.param_grid_low_memory_mode)
            self.assertEqual(params.param_grid_size, "small")
            self.assertEqual(params.parallel_fits, rules["combined_ml_pangraph"].params.parallel_fits)
            self.assertEqual(rules["combined_ml_mc"].script, rules["combined_ml"].script)
            self.assertEqual(rules["combined_ml_mc"].input.phenotype_table, rules["combined_ml"].input.phenotype_table)
            references = rules["mc_report"].params.reference_ml_dirs
            self.assertTrue(all("experiments/pangraph" in Path(p).as_posix() for p in references.values()))

    def test_experiment_rules_write_separately(self):
        with self.workflow() as (root, rules):
            output = root / "output"
            for name in ("mc_input", "mc_build", "mc_features", "merge_mc_set", "pivot_mc_set",
                         "combined_ml_mc", "mc_report"):
                with self.subTest(rule=name):
                    for path in rules[name].output:
                        self.assertTrue(
                            Path(path).is_relative_to(output / "minigraph_cactus")
                            or Path(path).is_relative_to(output / "experiments" / "minigraph_cactus"),
                            path,
                        )

    def run_input_rule(self, root, rules, reference_text):
        store = root / "store"
        store.mkdir()
        (store / SAMPLE).write_text(">contig_1 len=8\nacgtacgt\n>contig_2\nTTTT\n")
        (store / OTHER).write_text(">c\nGGGG\n")
        reference = root / "reference"
        reference.write_text(reference_text)
        command = format_command(
            rules["mc_input"].shellcmd,
            input=SimpleNamespace(samples=[str(store / OTHER), str(store / SAMPLE)], reference=str(reference)),
            output=SimpleNamespace(samples=str(root / "samples.txt"), seqfile=str(root / "in" / "seqfile.txt")),
            params=SimpleNamespace(ref_name="MG1655", input_dir=str(root / "in")),
            log=str(root / "log"),
            PREPARE_MC_INPUT=str(REPO / "snakefiles" / "scripts" / "prepare_mc_input.py"),
        )
        (root / "in").mkdir()
        path = os.pathsep.join([str(Path(sys.executable).parent), os.environ["PATH"]])
        return subprocess.run(["bash", "-c", command], env={**os.environ, "PATH": path})

    def test_input_rule_writes_the_seqfile_with_the_reference_first(self):
        with self.workflow() as (root, rules):
            self.assertEqual(rules["mc_input"].input.reference, rules["snippy_runner"].input.reference)
            genbank = (
                "LOCUS       NC_000913               12 bp    DNA     circular CON 09-MAR-2022\n"
                "ORIGIN\n"
                "        1 agcttttcat tc\n"
                "//\n"
            )
            self.assertEqual(self.run_input_rule(root, rules, genbank).returncode, 0)
            seqfile = [line.split("\t") for line in (root / "in" / "seqfile.txt").read_text().splitlines()]
            self.assertEqual([name for name, _ in seqfile], ["MG1655", SAMPLE, OTHER])
            self.assertEqual(Path(seqfile[0][1]).read_text(), ">NC_000913\nAGCTTTTCATTC\n")
            self.assertEqual(Path(seqfile[1][1]).read_text(), ">contig_1\nACGTACGT\n>contig_2\nTTTT\n")

    def test_input_rule_rejects_a_protein_reference(self):
        with self.workflow() as (root, rules):
            protein = ">tr|A0A024L2A9|A0A024L2A9_ECOLX Protein\nMTVQTSKNPQVDI\n"
            self.assertNotEqual(self.run_input_rule(root, rules, protein).returncode, 0)
            self.assertIn("non-nucleotide characters", (root / "log").read_text())


if __name__ == "__main__":
    unittest.main()

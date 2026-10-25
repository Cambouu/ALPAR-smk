"""Checks for the PanGraph experiment rules (pangraph.smk); no tools or jobs are run."""

from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import gzip
import os
import subprocess
import sys
import unittest

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ConfigSettings, ResourceSettings
from snakemake.utils import format as format_command


REPO = Path(__file__).resolve().parents[1]
SAMPLE = "0123456789abcdef0123456789abcdef01234567"


class PangraphWorkflowTest(unittest.TestCase):
    @contextmanager
    def workflow(self, **overrides):
        with TemporaryDirectory(prefix="alpar-pangraph-") as temp:
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

    def test_feature_sets_merge_the_baseline_snp_and_gpa_tables(self):
        with self.workflow() as (_, rules):
            merge = rules["merge_pangraph_set"].input[0]
            self.assertEqual(merge(SimpleNamespace(pangraph_set="snp_gpa")), [
                rules["binary_mutation_table"].output[0], rules["binary_gpa"].output[0],
            ])
            self.assertEqual(merge(SimpleNamespace(pangraph_set="gpa_pangraph")), [
                rules["binary_gpa"].output[0], rules["pangraph_features"].output.features,
            ])
            self.assertEqual(rules["pivot_pangraph_set"].shellcmd, rules["pivot_merged_features_miller"].shellcmd)

    def test_configured_sets_and_filter_reach_the_rules(self):
        with self.workflow(
            pangraph_feature_sets={"pangraph": ["pangraph"], "snp_gpa_pangraph": ["snp", "gpa", "pangraph"]},
            pangraph_min_carriers=3,
        ) as (_, rules):
            self.assertEqual(rules["pangraph_features"].params.min_carriers, 3)
            self.assertIn("--min-carriers {params.min_carriers}", rules["pangraph_features"].shellcmd)
            # Snakemake merges config dictionaries with config.yaml's, so the default sets stay.
            self.assertEqual(set(rules["pangraph_report"].params.feature_sets), {
                "snp_gpa", "gpa_pangraph", "pangraph", "snp_gpa_pangraph",
            })
            self.assertEqual(rules["pangraph_report"].params.reference_set, "snp_gpa")
            self.assertEqual(rules["pangraph_report"].params.feature_sets["snp_gpa_pangraph"], ["snp", "gpa", "pangraph"])
            results = rules["pangraph_report"].input.results
            self.assertEqual(len(results), 4)  # four sets, one antibiotic
            self.assertTrue(all("/drug/" in path.replace("\\", "/") for path in results))

    def test_ml_matches_the_baseline_search_for_every_set(self):
        with self.workflow() as (_, rules):
            params = rules["combined_ml_pangraph"].params
            self.assertTrue(params.param_grid_low_memory_mode)
            self.assertEqual(params.param_grid_size, "small")
            self.assertEqual(params.parallel_fits, rules["combined_ml"].params.parallel_fits)
            self.assertEqual(rules["combined_ml_pangraph"].script, rules["combined_ml"].script)
            self.assertEqual(
                rules["combined_ml_pangraph"].input.phenotype_table,
                rules["combined_ml"].input.phenotype_table,
            )

    def test_pangraph_build_uses_the_pinned_image_and_verifies(self):
        with self.workflow() as (_, rules):
            build = rules["pangraph_build"]
            self.assertIsNone(build.conda_env)
            self.assertTrue(build.container_img.startswith("docker://docker.io/cambouu/alpar-smk-pangraph@sha256:"))
            self.assertIn("pangraph build {input:q} -o {output:q} -j {threads} --verify", build.shellcmd)
            self.assertEqual(build.resources["_cores"], 8)

    def test_experiment_rules_write_separately(self):
        with self.workflow() as (root, rules):
            output = root / "output"
            names = {
                "pangraph_input", "pangraph_build", "pangraph_features", "merge_pangraph_set",
                "pivot_pangraph_set", "combined_ml_pangraph", "pangraph_report",
            }
            for name in names:
                with self.subTest(rule=name):
                    for path in rules[name].output:
                        self.assertTrue(
                            Path(path).is_relative_to(output / "pangraph")
                            or Path(path).is_relative_to(output / "experiments" / "pangraph"),
                            path,
                        )

    def run_input_rule(self, root, rules, reference_text):
        store = root / "store"
        store.mkdir()
        (store / SAMPLE).write_text(">contig_1 len=8\nacgtacgt\n>contig_2\nTTTT\n")
        reference = root / "reference"
        reference.write_text(reference_text)
        command = format_command(
            rules["pangraph_input"].shellcmd,
            input=SimpleNamespace(samples=[str(store / SAMPLE)], reference=str(reference)),
            output=SimpleNamespace(fasta=str(root / "input.fa.gz"), samples=str(root / "samples.txt")),
            params=SimpleNamespace(ref_name="MG1655"),
            log=str(root / "log"),
            PREPARE_PANGRAPH_INPUT=str(REPO / "snakefiles" / "scripts" / "prepare_pangraph_input.py"),
        )
        # The rule's image provides `python`; a login node may only have python3.
        path = os.pathsep.join([str(Path(sys.executable).parent), os.environ["PATH"]])
        return subprocess.run(["bash", "-c", command], env={**os.environ, "PATH": path})

    def test_input_rule_writes_pansn_names_for_every_sample(self):
        with self.workflow() as (root, rules):
            # The GBFF snippy uses, not fasta_file (Prokka's protein database).
            self.assertEqual(rules["pangraph_input"].input.reference, rules["snippy_runner"].input.reference)
            genbank = (
                "LOCUS       NC_000913               12 bp    DNA     circular CON 09-MAR-2022\n"
                "DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.\n"
                "FEATURES             Location/Qualifiers\n"
                "     source          1..12\n"
                "ORIGIN\n"
                "        1 agcttttcat tc\n"
                "//\n"
            )
            self.assertEqual(self.run_input_rule(root, rules, genbank).returncode, 0)
            self.assertEqual((root / "samples.txt").read_text().split(), [str(root / "store" / SAMPLE)])
            with gzip.open(root / "input.fa.gz", "rt") as handle:
                lines = [line.strip() for line in handle]
            self.assertEqual(lines, [
                ">MG1655#1#NC_000913", "AGCTTTTCATTC",
                f">{SAMPLE}#1#contig_1", "ACGTACGT", f">{SAMPLE}#1#contig_2", "TTTT",
            ])

    def test_input_rule_rejects_a_protein_reference(self):
        with self.workflow() as (root, rules):
            protein = ">tr|A0A024L2A9|A0A024L2A9_ECOLX Protein\nMTVQTSKNPQVDI\n"
            self.assertNotEqual(self.run_input_rule(root, rules, protein).returncode, 0)
            self.assertIn("non-nucleotide characters", (root / "log").read_text())


if __name__ == "__main__":
    unittest.main()

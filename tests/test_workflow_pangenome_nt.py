"""Checks for the nucleotide pangenome experiment rules (pangenome_nt.smk)."""

from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import subprocess
import unittest

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ConfigSettings, ResourceSettings


REPO = Path(__file__).resolve().parents[1]
SAMPLE = "0123456789abcdef0123456789abcdef01234567"


class NucleotidePangenomeWorkflowTest(unittest.TestCase):
    @contextmanager
    def workflow(self):
        with TemporaryDirectory(prefix="alpar-pangenome-nt-") as temp:
            root = Path(temp)
            (root / "input" / "drug").mkdir(parents=True)
            with SnakemakeApi() as api:
                workflow_api = api.workflow(
                    ResourceSettings(cores=2),
                    config_settings=ConfigSettings(config={
                        "input_dir": str(root / "input"),
                        "output_dir": str(root / "output"),
                        "env_dir": None,
                    }),
                    snakefile=REPO / "Snakefile",
                    workdir=root,
                )
                yield root, {rule.name: rule for rule in workflow_api._workflow.rules}

    def test_gene_names_match_the_clustered_protein_names(self):
        with self.workflow() as (root, rules):
            prokka = root / SAMPLE
            prokka.mkdir()
            (prokka / f"{SAMPLE}.faa").write_text(
                ">ABCDEFGH_00001 DNA gyrase subunit A (EC 5.6.2.2)\nMSDL\n"
                ">ABCDEFGH_00002 hypothetical protein\nMKT\n"
            )
            (prokka / f"{SAMPLE}.ffn").write_text(
                ">ABCDEFGH_00001 DNA gyrase subunit A (EC 5.6.2.2)\nATGAGCGACCTT\n"
                ">ABCDEFGH_00002 hypothetical protein\nATGAAAACC\n"
                ">ABCDEFGH_00003 tRNA-Leu(caa)\nGCCGAAG\n"
            )
            proteins, genes = root / "proteins.faa", root / "genes.ffn"
            for rule, kwargs in (
                ("cdhit_protein_name_corrector", dict(
                    input=prokka / f"{SAMPLE}.faa", output=proteins,
                    wildcards=SimpleNamespace(sample=SAMPLE),
                )),
                ("combine_ffn_files", dict(input=prokka / f"{SAMPLE}.ffn", output=genes)),
            ):
                subprocess.run(["bash", "-c", rules[rule].shellcmd.format(**kwargs)], check=True)

            def headers(path):
                return [line for line in path.read_text().splitlines() if line.startswith(">")]

            self.assertEqual(headers(proteins), [
                f">{SAMPLE}_ABCDEFGH_00001_DNA_gyrase_subunit_A__EC_5_6_2_2",
                f">{SAMPLE}_ABCDEFGH_00002_hypothetical_protein",
            ])
            # Every protein has its coding sequence; RNA genes are extra and unclustered.
            self.assertEqual(headers(genes)[:2], headers(proteins))
            self.assertEqual(len(headers(genes)), 3)
            self.assertIn("ATGAGCGACCTT", genes.read_text())

    def test_nucleotide_rules_reuse_gpa_clusters_and_write_separately(self):
        with self.workflow() as (root, rules):
            split = rules["split_cluster_fasta_nt"]
            self.assertEqual(split.input.cdhit_clstr, rules["cdhit_runner"].output.clstr)
            self.assertEqual(split.input.combined_proteins, rules["combine_ffn_files"].output[0])
            self.assertEqual(
                rules["split_cluster_fasta"].input.combined_proteins,
                rules["combine_faa_files"].output[0],
            )
            output = root / "output"
            for name, rule in rules.items():
                if name.endswith("_nt") and name not in {"pangenome_nt", "experiment_pangenome_nt"}:
                    with self.subTest(rule=name):
                        for path in rule.output:
                            self.assertTrue(
                                Path(path).is_relative_to(output / "pangenome_nt")
                                or Path(path).is_relative_to(output / "experiments" / "snp_gpa_pangenome_nt"),
                                path,
                            )

    def test_nucleotide_alignment_uses_the_dna_matrix_not_translation(self):
        with self.workflow() as (_, rules):
            rule = rules["panpa_align_nt"]
            self.assertEqual(rule.params.sub_matrix, "dna")
            self.assertIn('--sub_matrix "$SUB_MATRIX"', rule.shellcmd)
            self.assertNotIn("--dna", rule.shellcmd)
            self.assertNotIn("--sub_matrix", rules["panpa_align"].shellcmd)

    def test_only_nucleotide_bubble_features_cap_paths_and_logging(self):
        with self.workflow() as (_, rules):
            nucleotide = rules["bubble_features_nt"].shellcmd
            self.assertIn('--max-bubble-paths "$PATH_LIMIT" --log-level WARNING', nucleotide)
            self.assertIn("PATH_LIMIT={BUBBLE_PATH_LIMIT}", nucleotide)
            protein = rules["bubble_features"].shellcmd
            self.assertNotIn("--max-bubble-paths", protein)
            self.assertNotIn("--log-level", protein)


if __name__ == "__main__":
    unittest.main()

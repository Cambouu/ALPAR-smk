"""Checks for the CD-HIT variant experiment (cdhit_variants.smk)."""

import ast
from contextlib import contextmanager
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ConfigSettings, ResourceSettings


REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "snakefiles" / "scripts"
A = "a" * 40
B = "b" * 40
C = "c" * 40

# Member lines as cd-hit 4.8.1 writes them; -G 0 adds the alignment coordinates.
STANDARD_LINE = f"1\t290aa, >{B}_BBBBBBBB_00007_DNA_gyrase_subunit_A... at 95.00%"
LOCAL_IDENTITY_LINE = f"1\t6324aa, >{B}_DMLIMDBC_02828_Ig_like_domain_repeat_protein... at 72:6324:958:7227/84.48%"
REPRESENTATIVE_LINE = f"0\t300aa, >{A}_AAAAAAAA_00001_DNA_gyrase_subunit_A... *"


def module_constant(path, name):
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {path}")


GPA_REGEX = re.compile(module_constant(SCRIPTS / "binary_gpa_cdhit.py", "CLSTR_MEMBER_REGEX"))
SPLIT_REGEX = re.compile(module_constant(SCRIPTS / "split_cluster_fasta.py", "PROTEIN_RAW_REGEX"))


class ClstrParsingTest(unittest.TestCase):
    def test_both_parsers_accept_every_member_line_format(self):
        for line in (STANDARD_LINE, LOCAL_IDENTITY_LINE, REPRESENTATIVE_LINE):
            with self.subTest(line=line):
                gpa = GPA_REGEX.match(line)
                self.assertIsNotNone(gpa)
                self.assertEqual(len(gpa.group("strain")), 40)
                self.assertEqual(gpa.group("is_ref") is not None, line.endswith("*"))
                split = SPLIT_REGEX.match(line)
                self.assertIsNotNone(split)
                self.assertTrue(split.group("protein").startswith(">"))

    def test_merge_cdhit_steps_joins_first_clusters_under_second(self):
        first = (
            ">Cluster 0\n"
            f"0\t300aa, >{A}_AAAAAAAA_00001_DNA_gyrase_subunit_A... *\n"
            f"1\t290aa, >{B}_BBBBBBBB_00007_DNA_gyrase_subunit_A... at 95.00%\n"
            ">Cluster 1\n"
            f"0\t200aa, >{C}_CCCCCCCC_00002_hypothetical_protein... *\n"
            ">Cluster 2\n"
            f"0\t150aa, >{A}_AAAAAAAA_00003_putative_protein... *\n"
            f"1\t140aa, >{C}_CCCCCCCC_00004_putative_protein... at 1:140:3:142/72.00%\n"
        )
        second = (
            ">Cluster 0\n"
            f"0\t300aa, >{A}_AAAAAAAA_00001_DNA_gyrase_subunit_A... *\n"
            f"1\t150aa, >{A}_AAAAAAAA_00003_putative_protein... at 61.50%\n"
            ">Cluster 1\n"
            f"0\t200aa, >{C}_CCCCCCCC_00002_hypothetical_protein... *\n"
        )
        with TemporaryDirectory() as temp:
            temp = Path(temp)
            (temp / "first.clstr").write_text(first)
            (temp / "second.clstr").write_text(second)
            subprocess.run(
                [sys.executable, SCRIPTS / "merge_cdhit_steps.py",
                 temp / "first.clstr", temp / "second.clstr", temp / "out.clstr"],
                check=True, capture_output=True,
            )
            merged = (temp / "out.clstr").read_text().splitlines()
        self.assertEqual(merged, [
            ">Cluster 0",
            f"0\t300aa, >{A}_AAAAAAAA_00001_DNA_gyrase_subunit_A... *",
            f"1\t290aa, >{B}_BBBBBBBB_00007_DNA_gyrase_subunit_A... at 95.00%",
            f"2\t150aa, >{A}_AAAAAAAA_00003_putative_protein... at 61.50%",
            f"3\t140aa, >{C}_CCCCCCCC_00004_putative_protein... at 1:140:3:142/72.00%",
            ">Cluster 1",
            f"0\t200aa, >{C}_CCCCCCCC_00002_hypothetical_protein... *",
        ])
        for line in merged:
            if not line.startswith(">"):
                self.assertIsNotNone(GPA_REGEX.match(line), line)

    def test_merge_cdhit_steps_fails_when_a_first_cluster_is_missing(self):
        first = (
            ">Cluster 0\n"
            f"0\t300aa, >{A}_AAAAAAAA_00001_x... *\n"
            ">Cluster 1\n"
            f"0\t200aa, >{C}_CCCCCCCC_00002_y... *\n"
        )
        second = ">Cluster 0\n" + f"0\t300aa, >{A}_AAAAAAAA_00001_x... *\n"
        with TemporaryDirectory() as temp:
            temp = Path(temp)
            (temp / "first.clstr").write_text(first)
            (temp / "second.clstr").write_text(second)
            result = subprocess.run(
                [sys.executable, SCRIPTS / "merge_cdhit_steps.py",
                 temp / "first.clstr", temp / "second.clstr", temp / "out.clstr"],
                capture_output=True, text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing", result.stderr)


class CdhitVariantWorkflowTest(unittest.TestCase):
    @contextmanager
    def workflow(self):
        with TemporaryDirectory(prefix="alpar-cdhit-variants-") as temp:
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

    def test_variant_rules_write_separately(self):
        with self.workflow() as (root, rules):
            output = root / "output"
            names = [n for n in rules if n.endswith("_variant") or "cdhit_variant" in n or n == "merge_cdhit_steps"]
            self.assertIn("combined_ml_cdhit_variant", names)
            for name in names:
                with self.subTest(rule=name):
                    for path in rules[name].output:
                        self.assertTrue(
                            Path(path).is_relative_to(output / "cdhit_variants")
                            or Path(path).is_relative_to(output / "experiments" / "cdhit_variants"),
                            path,
                        )

    def test_variant_panproteome_uses_the_fixed_bubble_features_command(self):
        with self.workflow() as (_, rules):
            self.assertEqual(rules["bubble_features_variant"].shellcmd, rules["bubble_features_nt"].shellcmd)
            self.assertEqual(rules["panpa_align_variant"].shellcmd, rules["panpa_align"].shellcmd)
            self.assertEqual(rules["align_clusters_variant"].shellcmd, rules["align_clusters"].shellcmd)


if __name__ == "__main__":
    unittest.main()

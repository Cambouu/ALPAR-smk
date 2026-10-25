"""Tests for pangraph_features.py on a hand-made PanGraph JSON with known answers.

Reference REF has one contig 'chr' (26 bp) made of three blocks:
  b100 '+' at chr 1-10   consensus ACGTACGTAC
  b200 '-' at chr 11-20  consensus AAAACCCCGG  (chr 11-20 reads CCGGGGTTTT)
  b500 '+' at chr 21-26  consensus ACGT, REF carries an insertion TT before column 2 (chr reads ACTTGT)
Strain S1: b100 with G>T at column 2, columns 5-6 deleted, GG inserted before column 9;
           b200 with A>G at column 0 (T>C at chr 20 in reference orientation); b400 unchanged.
Strain S2: b300 and b400 (no reference copy; b400 has G>A at column 3); b500 without REF's TT.

Run: python -m unittest tests.test_pangraph_features
"""

import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "snakefiles", "scripts", "pangraph_features.py")


def edits(subs=(), dels=(), inss=()):
    return {"subs": [{"pos": p, "alt": a} for p, a in subs],
            "dels": [{"pos": p, "len": n} for p, n in dels],
            "inss": [{"pos": p, "seq": s} for p, s in inss]}


def toy_graph():
    nodes = {
        1: (100, 0, "+", [0, 10]), 2: (200, 0, "-", [10, 20]), 8: (500, 0, "+", [20, 26]),
        3: (100, 1, "+", [0, 10]), 4: (200, 1, "+", [10, 20]),
        5: (300, 2, "+", [0, 8]), 6: (400, 2, "+", [8, 16]), 9: (500, 2, "+", [16, 20]),
        7: (400, 3, "+", [0, 8]),
    }
    blocks = {
        100: ("ACGTACGTAC", {1: edits(), 3: edits(subs=[(2, "T")], dels=[(5, 2)], inss=[(9, "GG")])}),
        200: ("AAAACCCCGG", {2: edits(), 4: edits(subs=[(0, "G")])}),
        300: ("TTTTGGGG", {5: edits()}),
        400: ("CATGCATG", {6: edits(subs=[(3, "A")]), 7: edits()}),
        500: ("ACGT", {8: edits(inss=[(2, "TT")]), 9: edits()}),
    }
    paths = {0: ("REF#1#chr", [1, 2, 8], 26), 1: ("S1#1#c1", [3, 4], 20),
             2: ("S2#1#c1", [5, 6, 9], 20), 3: ("S1#1#c2", [7], 8)}
    return {
        "paths": {str(i): {"id": i, "name": n, "nodes": ns, "tot_len": L, "circular": False}
                  for i, (n, ns, L) in paths.items()},
        "blocks": {str(b): {"id": b, "consensus": c, "alignments": {str(k): v for k, v in al.items()}}
                   for b, (c, al) in blocks.items()},
        "nodes": {str(i): {"id": i, "block_id": b, "path_id": p, "strand": s, "position": pos}
                  for i, (b, p, s, pos) in nodes.items()},
    }


class PangraphFeaturesTest(unittest.TestCase):
    def run_extractor(self, *extra, carriers=None):
        with tempfile.TemporaryDirectory() as tmp:
            graph = os.path.join(tmp, "graph.json.gz")
            with gzip.open(graph, "wt") as handle:
                json.dump(toy_graph(), handle)
            out = os.path.join(tmp, "features.tsv")
            if carriers is not None:
                extra += ("--carriers", os.path.join(tmp, "carriers.tsv"))
            subprocess.run([sys.executable, SCRIPT, graph, out, "--ref", "REF", *extra],
                           check=True, capture_output=True, text=True)
            with open(out) as handle:
                rows = [line.rstrip("\n").split("\t") for line in handle]
            if carriers is not None:
                with open(os.path.join(tmp, "carriers.tsv")) as handle:
                    carriers.update(line.rstrip("\n").split("\t") for line in handle)
        features = {}
        for strain, feature, value in rows:
            self.assertEqual(value, "1")
            features.setdefault(strain, set()).add(feature)
        return features

    def test_expected_features(self):
        self.assertEqual(self.run_extractor(), {
            "S1": {
                "pangraph|blk|b100|chr:1-10",
                "pangraph|snp|3|G>T",
                "pangraph|del|6-7",
                "pangraph|ins|9^10|GG",
                "pangraph|blk|b200|chr:11-20",
                "pangraph|snp|20|T>C",
                "pangraph|blk|b400",
            },
            "S2": {
                "pangraph|blk|b300",
                "pangraph|blk|b400",
                "pangraph|snp|b400:3|G>A",
                "pangraph|blk|b500|chr:21-26",
                "pangraph|del|22^25|TT",
            },
        })

    def test_reference_not_written_and_rename(self):
        with tempfile.TemporaryDirectory() as tmp:
            mapping = os.path.join(tmp, "rename.tsv")
            with open(mapping, "w") as handle:
                handle.write("S1\tchecksum_one\n")
            features = self.run_extractor("--rename", mapping)
        self.assertEqual(set(features), {"checksum_one", "S2"})

    def test_reference_bases_match_reference_sequence(self):
        # chr = ACGTACGTAC + revcomp(AAAACCCCGG) + ACTTGT
        chr_seq = "ACGTACGTAC" + "CCGGGGTTTT" + "ACTTGT"
        for feature in set().union(*self.run_extractor().values()):
            _, kind, rest = feature.split("|", 2)
            if kind == "snp" and not rest.startswith("b"):
                pos, alleles = rest.split("|")
                self.assertEqual(chr_seq[int(pos) - 1], alleles.split(">")[0], feature)

    def test_carriers_count_strains_and_filter_keeps_shared_features(self):
        carriers = {}
        self.assertEqual(
            self.run_extractor("--min-carriers", "2", carriers=carriers),
            {"S1": {"pangraph|blk|b400"}, "S2": {"pangraph|blk|b400"}},
        )
        # Every feature is counted before the filter: S1's 7 and S2's 5 share only b400.
        self.assertEqual(len(carriers), 11)
        self.assertEqual(carriers["pangraph|blk|b400"], "2")
        self.assertEqual(set(carriers.values()), {"1", "2"})

    def test_default_min_carriers_keeps_every_feature(self):
        self.assertEqual(self.run_extractor("--min-carriers", "1"), self.run_extractor())


if __name__ == "__main__":
    unittest.main()

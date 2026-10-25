"""Per-strain features from a PanGraph (v1 JSON), in ALPAR's long format: strain<TAB>feature<TAB>1.

PanGraph splits the genomes into blocks (homologous segments with a consensus and each copy's
edits against it); each input contig is a path of blocks. Path names must be <strain>#<hap>#<contig>.

Feature families (all prefixed "pangraph|" so they cannot collide with other ALPAR families):

  pangraph|blk|b<block>[|<REF>:<start>-<end>]
      the strain has at least one copy of the block (accessory genes, plasmids, islands ...).
      The reference range is added when the reference genome has exactly one copy of the block.

  Variants of each strain copy against an anchor row:
    - the reference copy, if the block has exactly one copy on the reference: names use the
      reference's 1-based coordinates and bases in the reference orientation, like snippy:
        pangraph|snp|<pos>|<R>><A>          substitution
        pangraph|del|<start>-<end>          reference bases (aligned to the consensus) missing in the strain
        pangraph|ins|<left>^<right>|<seq>   strain bases between two reference positions
        pangraph|del|<left>^<right>|<seq>   reference bases that are an insertion vs the consensus, missing in the strain
        pangraph|cpx|<left>^<right>|<R>><A> such reference bases replaced by different inserted bases
    - otherwise the block consensus: the same types with b<block>:<column> instead of a
      reference position (consensus orientation, 0-based columns).
  Sequences longer than 50 bp are named by length and an md5 prefix instead of their sequence.

A strain with several copies of a block carries the union of its copies' variants. The reference
itself is used only as the anchor and is not written to the output.

Usage: python pangraph_features.py graph.json.gz out.tsv [--ref MG1655] [--rename map.tsv]
                                   [--min-carriers N] [--carriers carriers.tsv]
       map.tsv: strain<TAB>output_id (e.g. ALPAR checksum); strains not in it keep their name.
       --min-carriers: write only features carried by at least N strains (label-free filter;
                       the reference is not counted). Default 1 keeps every feature.
       --carriers: also write feature<TAB>carriers for every feature, before the filter.
"""

import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")
MAX_SEQ = 50


def strain_of(path_name):
    return path_name.split("#", 1)[0]


def edits(block, node_id):
    """Column -> allele ('-' for deletion) and column -> inserted sequence (before that column)."""
    a = block["alignments"][str(node_id)]
    allele = {s["pos"]: s["alt"] for s in a["subs"]}
    for d in a["dels"]:
        for c in range(d["pos"], d["pos"] + d["len"]):
            allele[c] = "-"
    ins = {i["pos"]: i["seq"] for i in a["inss"]}
    return allele, ins


def column_to_ref(block, node):
    """For the reference copy: consensus column -> 1-based reference coordinate (aligned columns only)."""
    allele, ins = edits(block, node["id"])
    start, end = node["position"]
    length = end - start
    col_pos, idx = {}, 0
    for c in range(len(block["consensus"]) + 1):
        idx += len(ins.get(c, ""))
        if c < len(block["consensus"]) and allele.get(c) != "-":
            col_pos[c] = start + idx + 1 if node["strand"] == "+" else end - idx
            idx += 1
    assert idx == length, f"node {node['id']}: mapped {idx} bases, node spans {length}"
    return col_pos


def seq_name(seq):
    if len(seq) <= MAX_SEQ:
        return seq
    return f"len{len(seq)}_{hashlib.md5(seq.encode()).hexdigest()[:8]}"


def runs(columns, bridge):
    """Group sorted columns into runs; two columns join if every column between them is in `bridge`."""
    out = []
    for c in sorted(columns):
        if out and all(x in bridge for x in range(out[-1][-1] + 1, c)):
            out[-1].append(c)
        else:
            out.append([c])
    return out


def variants(block, node_id, anchor_id, col_pos, rc, label):
    """Feature names for one strain copy against the anchor (reference copy or consensus)."""
    cons = block["consensus"]
    s_all, s_ins = edits(block, node_id)
    a_all, a_ins = edits(block, anchor_id) if anchor_id is not None else ({}, {})
    out = set()
    dels, inss = [], []
    for c in sorted(set(s_all) | set(a_all)):
        s, a = s_all.get(c, cons[c]), a_all.get(c, cons[c])
        if s == a:
            continue
        if s != "-" and a != "-":
            if col_pos is not None:
                r, v = (a, s) if not rc else (a.translate(COMP), s.translate(COMP))
                out.add(f"pangraph|snp|{col_pos[c]}|{r}>{v}")
            else:
                out.add(f"pangraph|snp|{label}:{c}|{a}>{s}")
        elif s == "-":
            dels.append(c)
        else:
            inss.append(c)

    # gaps in both rows are not reference bases, so they do not split a deletion or insertion run
    both_gap = {c for c in a_all if a_all[c] == "-"}
    for run in runs(dels, both_gap):
        if col_pos is not None:
            ps = [col_pos[c] for c in run if c in col_pos]
            out.add(f"pangraph|del|{min(ps)}-{max(ps)}")
        else:
            out.add(f"pangraph|del|{label}:{run[0]}-{run[-1]}")

    def flank(c, step):
        while 0 <= c < len(cons):
            if c in col_pos:
                return col_pos[c]
            c += step
        return None

    def between(col):
        """'<lo>^<hi>' reference positions around a gap before `col`; '<pos>^edge' at a block end."""
        ends = sorted(p for p in (flank(col - 1, -1), flank(col, +1)) if p is not None)
        return f"{ends[0]}^{ends[1]}" if len(ends) == 2 else f"{ends[0]}^edge"

    def add_ins(col, seq):
        if col_pos is not None:
            s = seq if not rc else seq.translate(COMP)[::-1]
            out.add(f"pangraph|ins|{between(col)}|{seq_name(s)}")
        else:
            out.add(f"pangraph|ins|{label}:{col}|{seq_name(seq)}")

    s_gap = {c for c in s_all if s_all[c] == "-"}
    for run in runs(inss, s_gap):
        add_ins(run[0], "".join(s_all.get(c, cons[c]) for c in range(run[0], run[-1] + 1) if s_all.get(c, cons[c]) != "-"))
    for col in set(s_ins) | set(a_ins):
        s_seq, a_seq = s_ins.get(col, ""), a_ins.get(col, "")
        if s_seq == a_seq:
            continue
        if s_seq and not a_seq:
            add_ins(col, s_seq)
            continue
        # the anchor has bases with no consensus column here: missing (or replaced) in the strain
        kind = "del" if not s_seq else "cpx"
        alleles = seq_name(a_seq) if not s_seq else f"{seq_name(a_seq)}>{seq_name(s_seq)}"
        if col_pos is None:
            out.add(f"pangraph|{kind}|{label}:ins{col}|{alleles}")
        else:
            if rc:
                alleles = seq_name(a_seq.translate(COMP)[::-1]) if not s_seq else \
                    f"{seq_name(a_seq.translate(COMP)[::-1])}>{seq_name(s_seq.translate(COMP)[::-1])}"
            out.add(f"pangraph|{kind}|{between(col)}|{alleles}")
    return out


def kind(feature):
    """blk, or <type>_ref / <type>_block for variants named by reference coordinate / block column."""
    _, t, where = feature.split("|", 2)
    if t == "blk":
        return "blk"
    return f"{t}_block" if where.startswith("b") else f"{t}_ref"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("graph")
    ap.add_argument("out")
    ap.add_argument("--ref", default="MG1655")
    ap.add_argument("--rename")
    ap.add_argument("--min-carriers", type=int, default=1)
    ap.add_argument("--carriers")
    args = ap.parse_args()

    rename = {}
    if args.rename:
        for line in open(args.rename):
            k, v = line.rstrip("\n").split("\t")[:2]
            rename[k] = v

    g = json.load(gzip.open(args.graph, "rt") if args.graph.endswith(".gz") else open(args.graph))
    path_strain = {int(k): strain_of(v["name"]) for k, v in g["paths"].items()}
    ref_contig = {int(k): v["name"].split("#")[-1] for k, v in g["paths"].items() if strain_of(v["name"]) == args.ref}
    by_block = defaultdict(list)
    for n in g["nodes"].values():
        by_block[n["block_id"]].append(n)

    features = defaultdict(set)  # strain -> feature names
    for bid, nodes in by_block.items():
        block = g["blocks"][str(bid)]
        label = f"b{bid}"
        ref_nodes = [n for n in nodes if path_strain[n["path_id"]] == args.ref]
        anchor, col_pos, rc, blk = None, None, False, f"pangraph|blk|{label}"
        if len(ref_nodes) == 1:
            anchor = ref_nodes[0]["id"]
            col_pos = column_to_ref(block, ref_nodes[0])
            rc = ref_nodes[0]["strand"] == "-"
            lo, hi = sorted((ref_nodes[0]["position"][0] + 1, ref_nodes[0]["position"][1]))
            blk += f"|{ref_contig[ref_nodes[0]['path_id']]}:{lo}-{hi}"
        for n in nodes:
            strain = path_strain[n["path_id"]]
            if strain == args.ref:
                continue
            features[strain].add(blk)
            features[strain] |= variants(block, n["id"], anchor, col_pos, rc, label)

    # carriers are counted over all strains, without labels, so the filter cannot leak phenotypes
    carriers = Counter(f for feats in features.values() for f in feats)
    if args.carriers:
        with open(args.carriers, "w") as out:
            for f in sorted(carriers):
                out.write(f"{f}\t{carriers[f]}\n")
    kept = {f for f, c in carriers.items() if c >= args.min_carriers}

    with open(args.out, "w") as out:
        for strain in sorted(features):
            sid = rename.get(strain, strain)
            for f in sorted(features[strain] & kept):
                out.write(f"{sid}\t{f}\t1\n")
    print(f"strains {len(features)}  distinct features {len(carriers)}  rows {sum(carriers.values())}")
    for k, v in sorted(Counter(kind(f) for f in carriers).items()):
        print(f"  {k:12s} {v}")
    if args.min_carriers > 1:
        print(f"min carriers {args.min_carriers}: kept features {len(kept)}  rows {sum(carriers[f] for f in kept)}")
        for k, v in sorted(Counter(kind(f) for f in kept).items()):
            print(f"  {k:12s} {v}")


if __name__ == "__main__":
    main()

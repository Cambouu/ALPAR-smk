"""Write a Minigraph-Cactus seqfile for a reference plus ALPAR assemblies.

The reference may be FASTA or GenBank (records named by their LOCUS name, as in snippy's ref.fa).

Every genome is rewritten as an uppercased FASTA (the PATRIC assemblies are lowercase, which cactus
treats as soft-masked) with contig names cut at the first whitespace. The sample name is the file
name (the ALPAR checksum); cactus reads '.' as a haplotype separator and '#' as a PanSN separator,
so names containing either are rejected. Assemblies are written in sorted order after the reference.

Usage: python prepare_mc_input.py --reference ref.fa --ref-name MG1655 --outdir DIR --seqfile out.txt
           ASSEMBLY... | @list.txt
"""

import argparse
import gzip
from pathlib import Path


# IUPAC nucleotide codes, checked after uppercasing: ALPAR's fasta_file is a protein database.
NUCLEOTIDES = frozenset("ACGTRYKMSWBDHVN")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def check_nucleotides(seq, src):
    if not NUCLEOTIDES.issuperset(seq):
        bad = "".join(sorted(set(seq) - NUCLEOTIDES))
        raise SystemExit(f"{src}: non-nucleotide characters {bad!r}")


def rewrite(src, dst):
    n_bases = 0
    with open_text(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                name = line[1:].split()[0]
                fout.write(f">{name}\n")
            else:
                seq = line.strip().upper()
                check_nucleotides(seq, src)
                n_bases += len(seq)
                fout.write(seq + "\n")
    return n_bases


def is_genbank(path):
    with open_text(path) as handle:
        return next((line for line in handle if line.strip()), "").startswith("LOCUS")


def write_genbank(src, dst):
    """GenBank (ALPAR's GBFF, snippy's reference) -> FASTA named by LOCUS, as in snippy's ref.fa."""
    n_bases, in_sequence = 0, False
    with open_text(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("LOCUS"):
                fout.write(f">{line.split()[1]}\n")
                in_sequence = False
            elif line.startswith("ORIGIN"):
                in_sequence = True
            elif line.startswith("//"):
                in_sequence = False
            elif in_sequence:
                seq = "".join(c for c in line if c.isalpha()).upper()
                check_nucleotides(seq, src)
                n_bases += len(seq)
                fout.write(seq + "\n")
    return n_bases


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--ref-name", default="MG1655")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seqfile", required=True)
    ap.add_argument("assemblies", nargs="+")
    a = ap.parse_args(argv)

    paths = []
    for item in a.assemblies:
        if item.startswith("@"):
            paths += [l.strip() for l in open(item[1:]) if l.strip()]
        else:
            paths.append(item)

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    samples = {}
    for p in paths:
        name = Path(p).name
        for suffix in (".gz", ".fna", ".fasta", ".fa"):
            name = name.removesuffix(suffix)
        if "." in name or "#" in name or name == a.ref_name:
            raise SystemExit(f"invalid sample name {name!r} for {p}")
        if name in samples:
            raise SystemExit(f"duplicate sample name {name!r}")
        samples[name] = p

    lines = []
    ref_out = outdir / f"{a.ref_name}.fa"
    if (write_genbank if is_genbank(a.reference) else rewrite)(a.reference, ref_out) == 0:
        raise SystemExit(f"empty reference {a.reference}")
    lines.append(f"{a.ref_name}\t{ref_out.resolve()}")
    for name in sorted(samples):
        dst = outdir / f"{name}.fa"
        if rewrite(samples[name], dst) == 0:
            raise SystemExit(f"empty assembly {samples[name]}")
        lines.append(f"{name}\t{dst.resolve()}")
    Path(a.seqfile).write_text("\n".join(lines) + "\n")
    print(f"wrote {len(lines)} genomes to {a.seqfile}")


if __name__ == "__main__":
    main()

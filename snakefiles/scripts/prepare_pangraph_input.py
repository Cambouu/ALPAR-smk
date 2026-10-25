"""Write the PanGraph input: one gzipped multi-FASTA with PanSN contig names.

Every contig becomes '<sample>#1#<contig>', where <sample> is the reference name for the
reference FASTA and the file name for each assembly (in ALPAR: the SHA-1 checksum under
data_preprocessing/data_checksum/). Sequences are uppercased: the local assemblies are entirely
lowercase, which some tools treat as soft-masked. Contig names are cut at the first whitespace.
PanGraph treats every record as its own path; the sample prefix is how the feature extractor
groups contigs into strains. The reference comes first, then the assemblies sorted by sample name
(PanGraph's output depends on input order). Standard library only (runs in ALPAR's python313 image).

The reference may be FASTA or GenBank (ALPAR's GBFF, the reference snippy uses; GenBank records are
named by their LOCUS name, as in snippy's ref.fa). Any sequence with a character outside PanGraph's
nucleotide alphabet stops the script: ALPAR's fasta_file is a protein database, not a genome.

Usage: python prepare_pangraph_input.py --reference ref.gbff --ref-name MG1655 --out input.fa.gz
                                        assembly [assembly ...] | @list.txt (one path per line)
"""

import argparse
import gzip
import os
import sys

# PanGraph's FASTA alphabet (IUPAC nucleotide codes), checked after uppercasing.
NUCLEOTIDES = frozenset("ACGTRYKMSWBDHVN")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def is_genbank(path):
    with opener(path) as handle:
        return next((line for line in handle if line.strip()), "").startswith("LOCUS")


def genbank_records(path):
    name, chunks, in_sequence = None, [], False
    with opener(path) as handle:
        for line in handle:
            if line.startswith("LOCUS"):
                name, chunks, in_sequence = line.split()[1], [], False
            elif line.startswith("ORIGIN"):
                in_sequence = True
            elif line.startswith("//"):
                yield name, "".join(chunks)
                in_sequence = False
            elif in_sequence:
                chunks.append("".join(c for c in line if c.isalpha()))


def records(path):
    if is_genbank(path):
        yield from genbank_records(path)
        return
    name, chunks = None, []
    with opener(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def sample_name(path):
    """File name without FASTA/compression extensions (ALPAR checksum files have none)."""
    name = os.path.basename(path)
    for ext in (".gz", ".fa", ".fna", ".fasta", ".fas"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    return name


def main():
    # @file reads one argument per line, for long assembly lists (e.g. 1,868 checksum paths)
    ap = argparse.ArgumentParser(fromfile_prefix_chars="@")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--ref-name", default="MG1655")
    ap.add_argument("--out", required=True)
    ap.add_argument("assemblies", nargs="+")
    args = ap.parse_args()

    # PanGraph's result depends on input order (two orders of the same 21 genomes shared 98.5% of
    # MG1655-coordinate SNP names), so sort by sample name: the same strain set always gives the same graph.
    samples = [(args.ref_name, args.reference)] + sorted((sample_name(p), p) for p in args.assemblies)
    seen = set()
    with gzip.open(args.out, "wt", compresslevel=1) as out:
        for sample, path in samples:
            if "#" in sample or sample in seen:
                sys.exit(f"invalid or duplicate sample name: {sample!r} ({path})")
            seen.add(sample)
            n = 0
            for contig, seq in records(path):
                if "#" in contig:
                    sys.exit(f"'#' in contig name {contig!r} ({path})")
                seq = seq.upper()
                if bad := set(seq) - NUCLEOTIDES:
                    sys.exit(f"non-nucleotide characters {''.join(sorted(bad))!r} in {contig!r} ({path})")
                out.write(f">{sample}#1#{contig}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
                n += 1
            if n == 0:
                sys.exit(f"no sequences in {path}")
    print(f"{len(samples)} samples written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

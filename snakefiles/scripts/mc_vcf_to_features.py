"""Minigraph-Cactus VCF (vg deconstruct, reference MG1655) -> ALPAR long feature table.

The VCF is split into biallelic records and atomised with bcftools (`norm -f REF -m -any -a`), so
complex alleles become SNPs and indels in reference coordinates, like Snippy's calls. A strain
carries a feature when any allele of its genotype is the ALT allele of the record. Missing genotypes
('.', the strain has no path through the site), the reference sample and cactus' own paths
(_MINIGRAPH_) give no feature.

Feature names (prefix `mc|`, so they cannot collide with other families):
  mc|snp|<pos>|<REF>><ALT>     single-base substitution
  mc|ind|<pos>|<REF>><ALT>     any other allele pair (indel, after atomisation)
Alleles longer than 50 bp are written as len<L>_<md5[:8]>. No name contains whitespace, '[', ']'
or '<', which XGBoost rejects.

Output: strain<TAB>feature<TAB>1, one line per strain and feature (sorted, C locale). Rows are
streamed to a temporary file and deduplicated with an external sort, so memory does not grow with
the number of strain-feature rows (about 10^8 for the whole collection).

Usage: python mc_vcf_to_features.py <in.vcf.gz> <ref.fa> <out.tsv> [--ref-sample MG1655]
"""

import argparse
import hashlib
import os
import subprocess
import sys


def allele_name(s):
    return s if len(s) <= 50 else f"len{len(s)}_{hashlib.md5(s.encode()).hexdigest()[:8]}"


def feature_name(pos, ref, alt):
    kind = "snp" if len(ref) == len(alt) == 1 else "ind"
    return f"mc|{kind}|{pos}|{allele_name(ref)}>{allele_name(alt)}"


def ref_contigs(ref_fa):
    return [l[1:].split()[0] for l in open(ref_fa) if l.startswith(">")]


def chrom_map(vcf, ref_fa):
    """Map VCF contig names such as MG1655#0#NC_000913 to the reference FASTA names."""
    header = subprocess.run(["bcftools", "view", "-h", vcf], capture_output=True, text=True, check=True).stdout
    vcf_contigs = [l.split("ID=")[1].split(",")[0].rstrip(">") for l in header.splitlines()
                   if l.startswith("##contig=<ID=")]
    names = set(ref_contigs(ref_fa))
    mapping = {}
    for c in vcf_contigs:
        if c in names:
            continue
        short = c.split("#")[-1]
        if short in names:
            mapping[c] = short
    return mapping


def samples_to_keep(vcf, ref_sample):
    """All sample columns except the reference and cactus' own paths such as _MINIGRAPH_."""
    names = subprocess.run(["bcftools", "query", "-l", vcf], capture_output=True, text=True,
                           check=True).stdout.split()
    return [s for s in names if s != ref_sample and not s.startswith("_")]


def records(vcf, ref_fa, keep_file):
    mapping = chrom_map(vcf, ref_fa)
    cmd = f"bcftools view -S {keep_file} -Ou {vcf}"
    if mapping:
        pairs = "".join(f"{k}\\t{v}\\n" for k, v in mapping.items())
        cmd += f" | bcftools annotate --rename-chrs <(printf '{pairs}') -Ou"
    cmd += (f" | bcftools norm -f {ref_fa} -m -any -a --atom-overlaps . -Ou"
            " | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n'")
    proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        yield line.rstrip("\n").split("\t")
    if proc.wait() != 0:
        raise SystemExit(f"bcftools pipeline failed: {cmd}")


def carriers(gts):
    for item in gts:
        sample, gt = item.rsplit("=", 1)
        if any(a == "1" for a in gt.replace("|", "/").split("/")):
            yield sample


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("ref_fa")
    ap.add_argument("out")
    ap.add_argument("--ref-sample", default="MG1655")
    ap.add_argument("--sort-memory", default="4G", help="buffer size of the final sort -u")
    a = ap.parse_args(argv)

    keep_file = a.out + ".samples"
    with open(keep_file, "w") as fh:
        fh.write("".join(f"{s}\n" for s in samples_to_keep(a.vcf, a.ref_sample)))

    unsorted = a.out + ".unsorted"
    distinct, strains, n_records = set(), set(), 0
    with open(unsorted, "w") as out:
        for f in records(a.vcf, a.ref_fa, keep_file):
            pos, ref, alt = f[0], f[1], f[2]
            if alt in ("*", ".") or ref == alt:
                continue
            n_records += 1
            name = feature_name(pos, ref, alt)
            for s in carriers(f[3:]):
                out.write(f"{s}\t{name}\t1\n")
                distinct.add(name)
                strains.add(s)
    tmpdir = os.environ.get("TMPDIR") or os.path.dirname(os.path.abspath(a.out))
    subprocess.run(["sort", "-u", "-S", a.sort_memory, "-T", tmpdir, "-o", a.out, unsorted],
                   check=True, env={**os.environ, "LC_ALL": "C"})
    os.remove(unsorted)
    with open(a.out) as fh:
        n_rows = sum(1 for _ in fh)
    print(f"{n_records} records, {len(distinct)} distinct features, {n_rows} rows, {len(strains)} strains",
          file=sys.stderr)


if __name__ == "__main__":
    main()

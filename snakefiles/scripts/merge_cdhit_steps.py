#!/usr/bin/env python3
"""Join a two-step CD-HIT clustering into one .clstr file.

Two-step (hierarchical) clustering from the CD-HIT user guide: FIRST clusters
every protein, SECOND clusters FIRST's representatives again at a lower
identity. Each SECOND cluster becomes one output cluster that holds every member
of the FIRST clusters whose representatives it contains, as cd-hit's
clstr_rev.pl does.

The output keeps SECOND's cluster numbers and representative (marked '*'). The
other FIRST representatives take their identity to it from SECOND; FIRST members
keep their identity to their FIRST representative. binary_gpa_cdhit.py and
split_cluster_fasta.py only read member names and the representative mark.

Usage: merge_cdhit_steps.py FIRST.clstr SECOND.clstr OUT.clstr
Standard library only. FIRST is indexed by byte offset, so memory stays small.
"""
import re
import sys

NAME = re.compile(rb">(\S+?)\.\.\. ")


def name_of(line):
    match = NAME.search(line)
    if match is None:
        raise SystemExit(f"no sequence name in {line!r}")
    return match.group(1)


def index_first(path):
    """Map each FIRST representative to its cluster's member lines as (offset, length)."""
    blocks = {}

    def add(rep, start, end):
        if rep is None:
            raise SystemExit(f"{path}: cluster at byte {start} has no representative")
        if rep in blocks:
            raise SystemExit(f"{path}: representative {rep!r} appears twice")
        blocks[rep] = (start, end - start)

    with open(path, "rb") as fh:
        offset = 0
        start = rep = None
        for line in fh:
            if line.startswith(b">Cluster"):
                if start is not None:
                    add(rep, start, offset)
                start, rep = offset + len(line), None
            elif line.rstrip().endswith(b"*"):
                rep = name_of(line)
            offset += len(line)
        if start is not None:
            add(rep, start, offset)
    return blocks


def second_clusters(path):
    """Yield (cluster number, member lines) for each SECOND cluster."""
    with open(path, "rb") as fh:
        number, members = None, []
        for line in fh:
            if line.startswith(b">Cluster"):
                if number is not None:
                    yield number, members
                number, members = int(line.split()[1]), []
            elif line.strip():
                members.append(line.rstrip(b"\r\n"))
        if number is not None:
            yield number, members


def main(first, second, out):
    blocks = index_first(first)
    n_first = len(blocks)
    n_clusters = n_members = 0
    with open(first, "rb") as src, open(out, "wb") as dst:
        for number, members in second_clusters(second):
            dst.write(b">Cluster %d\n" % number)
            n_clusters += 1
            k = 0
            for member in members:
                rep = name_of(member)
                is_second_rep = member.endswith(b"*")
                identity = None if is_second_rep else member.rsplit(b" at ", 1)[1]
                try:
                    start, length = blocks.pop(rep)
                except KeyError:
                    raise SystemExit(f"{rep!r} in {second} is not a representative in {first}, "
                                     "or appears twice in SECOND") from None
                src.seek(start)
                for line in src.read(length).splitlines():
                    if not line.strip():
                        continue
                    body = line.split(b"\t", 1)[1]
                    if body.endswith(b"*") and not is_second_rep:
                        body = body[:-1] + b"at " + identity
                    dst.write(b"%d\t%s\n" % (k, body))
                    k += 1
            n_members += k
    if blocks:
        raise SystemExit(f"{len(blocks)} of {n_first} FIRST clusters are missing from {second}, "
                         f"e.g. {next(iter(blocks))!r}")
    print(f"first_clusters={n_first} merged_clusters={n_clusters} members={n_members}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])

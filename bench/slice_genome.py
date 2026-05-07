#!/usr/bin/env python3
"""Slice a contiguous range out of a chromosome FASTA.

Reads bench/data_genomes/<asm>.<chr>.fa, takes a start position and a length,
writes a new FASTA file with the chosen subrange. Soft-masking and N-blocks
are preserved (we just slice the concatenated sequence string).

Used to make small workloads for fast profiler iteration:

    # First 10 Mbp of galGal6.chrZ + first 10 Mbp of taeGut2.chrZ
    python3 bench/slice_genome.py bench/data_genomes/galGal6.chrZ.fa  \\
                                  bench/data_genomes/galGal6.chrZ_0_10mb.fa \\
                                  --start 0 --length 10_000_000

We slice and rename the FASTA header so lastz sees a unique sequence id
(otherwise multiple slices of the same chromosome aliased to "chrZ" inside
lastz get confusing in MAF output).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> tuple[str, str]:
    """Returns (header_line_without_>, sequence_string).

    Assumes a single-sequence FASTA. We strip newlines from the body.
    """
    name: str | None = None
    parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    raise ValueError(f"{path}: multi-sequence FASTA not supported")
                name = line[1:].split(None, 1)[0]
            else:
                parts.append(line)
    if name is None:
        raise ValueError(f"{path}: no FASTA header")
    return name, "".join(parts)


def write_fasta(path: Path, name: str, seq: str, line_width: int = 50) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i : i + line_width])
            f.write("\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_fa", type=Path, help="Input FASTA (single sequence).")
    p.add_argument("output_fa", type=Path, help="Output FASTA path.")
    p.add_argument("--start", type=int, default=0,
                   help="0-based start position (default: 0).")
    p.add_argument("--length", type=int, required=True,
                   help="Number of bases to extract.")
    p.add_argument("--rename", type=str, default=None,
                   help="Rename the output sequence id (default: '<orig>_<start>_<end>').")
    args = p.parse_args(argv)

    name, seq = read_fasta(args.input_fa)
    n = len(seq)
    end = args.start + args.length
    if args.start < 0 or end > n:
        raise SystemExit(
            f"slice [{args.start}, {end}) out of bounds for {name} of length {n:,}")

    sub = seq[args.start : end]
    new_name = args.rename or f"{name}_{args.start}_{end}"

    args.output_fa.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.output_fa, new_name, sub)

    # Quick stats — useful in logs.
    upper = sum(1 for c in sub if c.isupper() and c != "N")
    lower = sum(1 for c in sub if c.islower())
    n_count = sum(1 for c in sub if c in ("N", "n"))
    print(f"wrote {args.output_fa}")
    print(f"  source:  {args.input_fa}  ({name}, {n:,} bp)")
    print(f"  slice:   [{args.start:,}, {end:,})  ({len(sub):,} bp)")
    print(f"  ACGT:    {upper:,}  ({100*upper/len(sub):.1f}%)")
    print(f"  acgt:    {lower:,}  ({100*lower/len(sub):.1f}%)  -- soft-masked repeats")
    print(f"  N/n:     {n_count:,}  ({100*n_count/len(sub):.1f}%)  -- assembly gaps")
    print(f"  header:  >{new_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate synthetic FASTA pairs for scalable lastz benchmarking.

We build a "target" sequence of random DNA, then derive a "query" from it by
applying a controlled mix of substitutions and short indels. The resulting
pair is similar enough that lastz finds non-trivial alignments (so we exercise
seed-and-extend, gapped extension, and chaining), but the divergence and size
are knobs we control.

Determinism: every output depends only on (length, divergence, indel_rate,
seed, name), so re-running with the same args yields byte-identical files.
That's important for benchmark reproducibility.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import sys
from pathlib import Path

ALPHABET = "ACGT"


def _rng(seed: int) -> random.Random:
    return random.Random(seed)


def random_dna(length: int, rng: random.Random) -> list[str]:
    # list[str] of length `length`, mutable so we can splice indels in.
    # Using random.choices is much faster than per-base random.choice.
    return rng.choices(ALPHABET, k=length)


def mutate(
    seq: list[str],
    rng: random.Random,
    sub_rate: float,
    indel_rate: float,
    max_indel: int = 5,
) -> list[str]:
    """Apply substitutions and short indels.

    sub_rate    : per-base probability of substitution
    indel_rate  : per-base probability of starting an indel event (length 1..max_indel)
    """
    out: list[str] = []
    i = 0
    n = len(seq)
    while i < n:
        r = rng.random()
        if r < sub_rate:
            base = seq[i]
            choices = [b for b in ALPHABET if b != base]
            out.append(rng.choice(choices))
            i += 1
        elif r < sub_rate + indel_rate:
            # 50/50 insertion vs deletion
            length = rng.randint(1, max_indel)
            if rng.random() < 0.5:
                out.extend(rng.choices(ALPHABET, k=length))
            else:
                i += length
        else:
            out.append(seq[i])
            i += 1
    return out


def write_fasta(path: Path, name: str, seq: list[str], width: int = 80) -> None:
    s = "".join(seq)
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for k in range(0, len(s), width):
            f.write(s[k : k + width])
            f.write("\n")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def make_pair(
    out_dir: Path,
    base_name: str,
    target_len: int,
    query_len: int | None,
    sub_rate: float,
    indel_rate: float,
    seed: int,
    force: bool = False,
) -> tuple[Path, Path, dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    target_path = out_dir / f"{base_name}.target.fa"
    query_path = out_dir / f"{base_name}.query.fa"

    meta = {
        "base_name": base_name,
        "target_len": target_len,
        "query_len_requested": query_len,
        "sub_rate": sub_rate,
        "indel_rate": indel_rate,
        "seed": seed,
    }

    if target_path.exists() and query_path.exists() and not force:
        meta["target_sha256"] = sha256_file(target_path)
        meta["query_sha256"] = sha256_file(query_path)
        meta["cached"] = True
        return target_path, query_path, meta

    rng_t = _rng(seed)
    target = random_dna(target_len, rng_t)

    rng_q = _rng(seed + 1)
    if query_len is None or query_len == target_len:
        query = mutate(target, rng_q, sub_rate, indel_rate)
    else:
        # Build a query of approximately query_len by tiling/truncating mutated copies.
        # Useful for read-mapping-ish scenarios where query >> or << target.
        query = []
        while len(query) < query_len:
            chunk = mutate(target, rng_q, sub_rate, indel_rate)
            query.extend(chunk)
        query = query[:query_len]

    write_fasta(target_path, f"{base_name}_target", target)
    write_fasta(query_path, f"{base_name}_query", query)

    meta["target_sha256"] = sha256_file(target_path)
    meta["query_sha256"] = sha256_file(query_path)
    meta["cached"] = False
    return target_path, query_path, meta


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=Path(__file__).parent / "data")
    p.add_argument("--name", required=True, help="Base name for output files")
    p.add_argument("--target-len", type=int, required=True)
    p.add_argument("--query-len", type=int, default=None,
                   help="Defaults to target-len (1:1 mutated copy)")
    p.add_argument("--sub-rate", type=float, default=0.05)
    p.add_argument("--indel-rate", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=0xC0FFEE)
    p.add_argument("--force", action="store_true", help="Regenerate even if files exist")
    args = p.parse_args(argv)

    t, q, meta = make_pair(
        out_dir=args.out_dir,
        base_name=args.name,
        target_len=args.target_len,
        query_len=args.query_len,
        sub_rate=args.sub_rate,
        indel_rate=args.indel_rate,
        seed=args.seed,
        force=args.force,
    )
    print(f"target: {t}  ({os.path.getsize(t):,} bytes, sha256={meta['target_sha256'][:12]}...)")
    print(f"query : {q}  ({os.path.getsize(q):,} bytes, sha256={meta['query_sha256'][:12]}...)")
    print(f"cached: {meta['cached']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

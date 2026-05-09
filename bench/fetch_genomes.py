#!/usr/bin/env python3
"""Fetch reference genome assemblies for the lastz benchmark harness.

Downloads UCSC 2bit files for the four assemblies SegAlign uses
(hg38, mm10, galGal6, taeGut2) into a scratch directory, verifies
md5 checksums against UCSC's published md5sum.txt, then creates
symlinks back into the workspace so workloads can find them.

Optionally extracts one chromosome per genome into a `slices/` directory.
We use a small in-tree 2bit reader (bench/twobit.py) instead of UCSC's
`twoBitToFa` because the prebuilt `twoBitToFa` binary requires glibc 2.33+
which isn't available on this host.

Layout produced:

    /scratch2/shiv1/lastz-bench-data/
        2bit/
            hg38.2bit, mm10.2bit, galGal6.2bit, taeGut2.2bit
            md5sum.<asm>.txt           # raw UCSC checksum manifests
        slices/
            hg38.chr19.fa              # ~58 Mbp, gene-rich human autosome
            mm10.chr10.fa              # ~130 Mbp, syntenic with hg38.chr19
            galGal6.chrZ.fa            # ~82 Mbp, bird sex chromosome
            taeGut2.chrZ.fa            # ~73 Mbp, paired with galGal6.chrZ

    bench/data_genomes/                # symlinks back into the workspace
        hg38.2bit -> /scratch2/.../2bit/hg38.2bit
        ... and so on for 2bits and slices ...

Idempotency: re-runs are no-ops if files exist and md5 matches. Pass --force
to redownload everything.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE_LINKS = HERE / "data_genomes"
SCRATCH_ROOT = Path("/scratch2/shiv1/lastz-bench-data")
UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath"

sys.path.insert(0, str(HERE))
from twobit import TwoBitReader  # noqa: E402


@dataclass
class Genome:
    asm: str               # short id (e.g. "hg38")
    species: str           # display name
    size_gbp: float        # informational
    slice_chr: str         # chromosome name to extract into slices/ as a sample
    notes: str = ""


# Same set SegAlign uses, so our wall-times can be directly compared to theirs.
# We also pull T2T-CHM13 (UCSC code "hs1") so we can run within-species
# (human-vs-human, ~99.9% identity) seed-weight experiments — the regime
# where bigger seeds + step>1 actually pay off.
GENOMES: list[Genome] = [
    Genome("hg38",    "human (Homo sapiens)",        3.08, "chr19",
           notes="GRCh38, Dec 2013. Gene-rich chr19 is mostly syntenic with mouse chr10."),
    Genome("hs1",     "human T2T-CHM13v2.0",         3.05, "chr1",
           notes="Jan 2022. Complete telomere-to-telomere human assembly. "
                 "Pairs with hg38 for human-vs-human (~99.9%% identity) experiments."),
    Genome("mm10",    "mouse (Mus musculus)",         2.72, "chr10",
           notes="GRCm38, Jan 2012. Matches what SegAlign benchmarks against."),
    Genome("rheMac10","rhesus macaque (Macaca mulatta)", 2.94, "chr19",
           notes="Mmul_10, Feb 2019. ~93% identity with hg38 in syntenic regions; "
                 "primate sister taxon, useful as a true mid-identity regime point "
                 "between hg-vs-CHM13 (~99.5%) and bird-Z (~75%)."),
    Genome("rn6",     "rat (Rattus norvegicus)",     2.87, "chr20",
           notes="Rnor_6.0, Jul 2014. ~85% identity with mm10 in syntenic regions; "
                 "rn6.chr20 is the syntenic counterpart of mm10.chr10. Useful as a "
                 "dense-synteny / mid-identity regime point to pin down the "
                 "seed->gapped bottleneck crossover (sits between bird-Z's 0.6 Mbp "
                 "aligned and rhesus's 23-37 Mbp aligned per 100 Mbp^2 product)."),
    Genome("galGal6", "chicken (Gallus gallus)",      1.05, "chrZ",
           notes="Mar 2018. Bird sex chr Z; birds are ZW (males ZZ, females ZW)."),
    Genome("taeGut2", "zebra finch (Taeniopygia guttata)", 1.02, "chrZ",
           notes="Feb 2013. Pairs with galGal6.chrZ for bird-vs-bird Z alignment."),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    """Run a shell command, raising on non-zero exit. Echoes the command."""
    print(f"  $ {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, check=True, **kw)


def md5_file(path: Path, block: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_md5sum_file(text: str) -> dict[str, str]:
    """UCSC md5sum.txt format: '<md5>  <filename>' per line."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out[parts[1]] = parts[0]
    return out


def fetch_text(url: str) -> str:
    r = subprocess.run(["curl", "-sS", "--max-time", "30", url],
                       check=True, capture_output=True, text=True)
    return r.stdout


def fetch_url_to_file(url: str, dest: Path, show_progress: bool = True) -> None:
    """Download `url` to `dest` atomically (write to .part, rename on success)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    args = ["curl", "-fL", "--retry", "3", "--retry-delay", "5",
            "-o", str(tmp), url]
    if show_progress:
        args.insert(1, "-#")  # progress bar to stderr
    else:
        args.insert(1, "-sS")
    sh(args)
    tmp.replace(dest)


def make_symlink(src: Path, dst: Path) -> None:
    """Create or refresh a symlink dst -> src."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_symlink() and Path(os.readlink(dst)) == src:
            return
        dst.unlink()
    dst.symlink_to(src)


# ---------------------------------------------------------------------------
# Fetch steps
# ---------------------------------------------------------------------------

def fetch_2bit(g: Genome, twobit_dir: Path, force: bool) -> tuple[Path, str]:
    """Download <asm>.2bit and verify against UCSC's md5sum.txt.

    Returns (path, md5).
    """
    out_path = twobit_dir / f"{g.asm}.2bit"
    md5_path = twobit_dir / f"md5sum.{g.asm}.txt"
    md5_url = f"{UCSC_BASE}/{g.asm}/bigZips/md5sum.txt"
    bit_url = f"{UCSC_BASE}/{g.asm}/bigZips/{g.asm}.2bit"

    # 1) Always re-fetch the (small) md5sum manifest so we know the expected hash.
    md5_text = fetch_text(md5_url)
    md5_path.write_text(md5_text)
    expected = parse_md5sum_file(md5_text).get(f"{g.asm}.2bit")
    if expected is None:
        raise RuntimeError(f"md5sum.txt for {g.asm} doesn't list {g.asm}.2bit")

    # 2) Skip download if the local file is already correct.
    if out_path.exists() and not force:
        actual = md5_file(out_path)
        if actual == expected:
            print(f"  {g.asm}.2bit: cached + md5 OK ({actual})")
            return out_path, actual
        print(f"  {g.asm}.2bit: md5 mismatch (have {actual}, want {expected}); refetching")

    # 3) Download.
    print(f"  {g.asm}.2bit: fetching {bit_url}")
    fetch_url_to_file(bit_url, out_path)

    # 4) Verify.
    actual = md5_file(out_path)
    if actual != expected:
        raise RuntimeError(
            f"md5 mismatch after download: have {actual}, want {expected} for {g.asm}.2bit"
        )
    print(f"  {g.asm}.2bit: download OK, md5 verified ({actual})")
    return out_path, actual


def extract_slice(g: Genome, twobit_path: Path, slices_dir: Path,
                  force: bool) -> Path:
    out = slices_dir / f"{g.asm}.{g.slice_chr}.fa"
    if out.exists() and not force:
        print(f"  slice {g.asm}.{g.slice_chr}: cached ({out.stat().st_size:,} B)")
        return out
    with TwoBitReader(twobit_path) as t:
        if not t.has(g.slice_chr):
            available = ", ".join(t.sequences()[:10])
            raise RuntimeError(
                f"chromosome {g.slice_chr} not found in {twobit_path}; "
                f"first 10 sequences are: {available}")
        n = t.write_fasta(g.slice_chr, out)
    print(f"  slice {g.asm}.{g.slice_chr}: extracted {n:,} bp "
          f"({out.stat().st_size:,} B)")
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scratch", type=Path, default=SCRATCH_ROOT,
                   help=f"Root for downloaded data (default: {SCRATCH_ROOT}).")
    p.add_argument("--links",   type=Path, default=WORKSPACE_LINKS,
                   help=f"Where to put symlinks back into the workspace "
                        f"(default: {WORKSPACE_LINKS}).")
    p.add_argument("--force",   action="store_true",
                   help="Re-download / re-extract everything.")
    p.add_argument("--no-slices", action="store_true",
                   help="Skip the per-genome slice extraction (only fetch 2bit).")
    p.add_argument("--only", nargs="+", choices=[g.asm for g in GENOMES],
                   help="Restrict to a subset of genomes.")
    args = p.parse_args(argv)

    scratch    = args.scratch
    twobit_dir = scratch / "2bit"
    slices_dir = scratch / "slices"
    for d in (twobit_dir, slices_dir):
        d.mkdir(parents=True, exist_ok=True)

    selected = [g for g in GENOMES if not args.only or g.asm in args.only]

    print("# 2bit downloads")
    fetched: list[tuple[Genome, Path, str]] = []
    for g in selected:
        path, md5 = fetch_2bit(g, twobit_dir, force=args.force)
        fetched.append((g, path, md5))

    extracted: list[tuple[Genome, Path]] = []
    if not args.no_slices:
        print()
        print("# slice extraction (via in-tree 2bit reader)")
        for g, twobit_path, _ in fetched:
            slice_path = extract_slice(g, twobit_path, slices_dir,
                                       force=args.force)
            extracted.append((g, slice_path))

    # Create / refresh symlinks back into the workspace.
    print()
    print(f"# symlinks → {args.links}")
    args.links.mkdir(parents=True, exist_ok=True)
    for g, twobit_path, _ in fetched:
        link = args.links / f"{g.asm}.2bit"
        make_symlink(twobit_path, link)
        print(f"  {link} -> {twobit_path}")
    for g, slice_path in extracted:
        link = args.links / slice_path.name
        make_symlink(slice_path, link)
        print(f"  {link} -> {slice_path}")

    print()
    print("# summary")
    for g, twobit_path, md5 in fetched:
        size_mb = twobit_path.stat().st_size / 1024 / 1024
        print(f"  {g.asm:8s} {g.species:35s} {size_mb:7.1f} MB  md5={md5}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

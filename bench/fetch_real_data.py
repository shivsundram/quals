#!/usr/bin/env python3
"""Fetch real biological sequences for the lastz benchmark harness.

Pulls a small, curated panel of homologous gene sequences via NCBI's
E-utilities (efetch) and writes them as one FASTA file per species into
bench/data_real/. No biopython dependency — just stdlib urllib.

Network policy:
- Files are cached on disk; re-runs without --force are no-ops.
- If a primary accession returns nothing, we try the listed fallbacks before
  giving up. This protects us from RefSeq accession bumps (predicted
  XM_* accessions in particular get versioned aggressively).
- We sleep 0.4 s between requests to stay within NCBI's anonymous rate cap
  (~3 req/sec).
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


@dataclass
class Target:
    species: str         # short id used in filenames (human, chimp, fly, ...)
    common: str          # display name
    latin: str           # Latin binomial
    gene: str            # gene symbol
    accessions: list[str]  # primary first, fallbacks after


# Curated panels. Each panel becomes a directory full of fa files.
PANELS: dict[str, list[Target]] = {
    # β-actin mRNA across human / chimp / fruit fly. ACTB in vertebrates is
    # ~1.8 kb and has 1:1 orthologs in human↔chimp; in Drosophila there are
    # six actin paralogs and Act5C is the canonical cytoplasmic one.
    "actb": [
        Target(
            species="human", common="human", latin="Homo sapiens",
            gene="ACTB",
            accessions=["NM_001101.5"],
        ),
        Target(
            species="chimp", common="chimp", latin="Pan troglodytes",
            gene="ACTB",
            # Curated RefSeq mRNA (looked up via NCBI esearch).
            accessions=[
                "NM_001009945.2",
                "NM_001009945.1",
            ],
        ),
        Target(
            species="fly", common="fruit fly", latin="Drosophila melanogaster",
            gene="Act5C",
            # Drosophila Act5C has five RefSeq transcript variants. Variant D
            # (1928 bp) is the closest in size to human/chimp ACTB (~1.8 kb),
            # which keeps the matrix cells size-comparable. Fall back to the
            # other variants if the version number bumps.
            accessions=[
                "NM_001014725.2",
                "NM_001014725.1",
                "NM_167053.2",  # variant A, 1637 bp
                "NM_078497.4",  # variant B, 2653 bp
            ],
        ),
    ],
}


def fetch_fasta(accession: str, retries: int = 3, delay: float = 1.0) -> str:
    qs = urlencode({"db": "nuccore", "id": accession,
                    "rettype": "fasta", "retmode": "text"})
    url = f"{EFETCH}?{qs}"
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            with urlopen(url, timeout=30) as r:
                body = r.read().decode()
            # NCBI returns an empty body or an error message (not HTTP 4xx)
            # when the accession is gone. Treat blank or non-FASTA as failure.
            if not body.strip().startswith(">"):
                raise RuntimeError(f"non-FASTA response from NCBI: {body[:120]!r}")
            return body
        except (HTTPError, URLError, RuntimeError) as e:
            last_err = e
            time.sleep(delay * (attempt + 1))
    raise RuntimeError(f"failed to fetch {accession} after {retries} attempts: {last_err}")


def fetch_target(target: Target) -> tuple[str, str]:
    """Try each accession in order; return (accession_used, fasta_text)."""
    last_err: Exception | None = None
    for acc in target.accessions:
        try:
            return acc, fetch_fasta(acc)
        except Exception as e:
            last_err = e
            print(f"    {acc}: {e}", file=sys.stderr)
    raise RuntimeError(
        f"all accessions failed for {target.species} {target.gene}: {last_err}"
    )


def normalize_fasta(text: str, species: str, gene: str, accession: str,
                    latin: str) -> str:
    """Rewrite the header to a uniform `>species_gene_accession latin name` form.

    This makes lastz alignment output much easier to read. The body is
    unchanged.
    """
    lines = text.split("\n", 1)
    body = lines[1] if len(lines) == 2 else ""
    new_header = f">{species}_{gene}_{accession} {latin}"
    return new_header + "\n" + body


def seq_length_from_fasta(text: str) -> int:
    # Skip every header line (starts with '>'); count nucleotide characters
    # in the body, matching the conventions of `wc` on a sequence-only FASTA.
    n = 0
    for line in text.splitlines():
        if line.startswith(">") or not line:
            continue
        n += sum(1 for c in line if c.upper() in "ACGTNRYSWKMBDHV")
    return n


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def fetch_panel(panel_name: str, out_dir: Path, force: bool, sleep_s: float) -> dict:
    if panel_name not in PANELS:
        raise SystemExit(f"unknown panel: {panel_name}; "
                         f"choices: {list(PANELS)}")
    panel = PANELS[panel_name]
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict = {"panel": panel_name, "entries": []}
    for t in panel:
        out_path = out_dir / f"{panel_name}_{t.species}.fa"
        if out_path.exists() and not force:
            text = out_path.read_text()
            entry = {
                "species": t.species, "gene": t.gene, "latin": t.latin,
                "path": str(out_path),
                "length_bp": seq_length_from_fasta(text),
                "sha256": sha256_text(text),
                "cached": True,
            }
            manifest["entries"].append(entry)
            print(f"  {t.species:6s} {t.gene:8s} cached: "
                  f"{entry['length_bp']:6d} bp → {out_path}")
            continue

        print(f"  {t.species:6s} {t.gene:8s} fetching ...", flush=True)
        used_acc, text = fetch_target(t)
        text = normalize_fasta(text, t.species, t.gene, used_acc, t.latin)
        out_path.write_text(text)
        length = seq_length_from_fasta(text)
        entry = {
            "species": t.species, "gene": t.gene, "latin": t.latin,
            "accession": used_acc,
            "path": str(out_path),
            "length_bp": length,
            "sha256": sha256_text(text),
            "cached": False,
        }
        manifest["entries"].append(entry)
        print(f"  {t.species:6s} {t.gene:8s} {used_acc:20s} "
              f"{length:6d} bp → {out_path}")
        time.sleep(sleep_s)
    return manifest


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel", default="actb", choices=PANELS.keys())
    p.add_argument("--out-dir", type=Path,
                   default=Path(__file__).parent / "data_real")
    p.add_argument("--force", action="store_true",
                   help="Re-fetch even if files are already cached.")
    p.add_argument("--sleep", type=float, default=0.4,
                   help="Seconds between NCBI requests (rate limit).")
    args = p.parse_args(argv)

    print(f"# panel: {args.panel}")
    print(f"# out:   {args.out_dir}")
    fetch_panel(args.panel, args.out_dir, args.force, args.sleep)
    return 0


if __name__ == "__main__":
    sys.exit(main())

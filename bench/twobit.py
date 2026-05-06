"""Minimal reader for UCSC .2bit genome assembly files.

Format spec: https://genome.ucsc.edu/FAQ/FAQformat.html#format7

Layout:
    Header (16 B):
        signature        u32  0x1A412743 (LE) or 0x4327411A (BE)
        version          u32  0
        sequenceCount    u32
        reserved         u32  0
    Index (sequenceCount entries):
        nameSize         u8
        name             nameSize ASCII bytes
        offset           u32  byte offset to this sequence's record
    Per-sequence record (at `offset`):
        dnaSize          u32  number of bases
        nBlockCount      u32
        nBlockStarts     u32 * nBlockCount
        nBlockSizes      u32 * nBlockCount
        maskBlockCount   u32
        maskBlockStarts  u32 * maskBlockCount
        maskBlockSizes   u32 * maskBlockCount
        reserved         u32  0
        packedDna        ceil(dnaSize/4) bytes; 2-bit codes: T=00, C=01, A=10, G=11
                         high bits first within each byte.

We expose just what the harness needs: list sequences, get a single sequence
as bytes, and write out one chromosome as FASTA. Soft-masking (lowercase for
maskBlock ranges) and N-blocks are honored.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import IO, NamedTuple

# 2-bit codes: T=00, C=01, A=10, G=11.
_CODES = b"TCAG"


class SeqIndexEntry(NamedTuple):
    name: str
    offset: int


class TwoBitReader:
    def __init__(self, path: Path):
        self.path = Path(path)
        self._fh: IO[bytes] = open(self.path, "rb")
        self._endian = "<"  # default little-endian; flipped below if needed
        self._read_header()
        self._index: dict[str, int] = {}
        self._read_index()

    # --- header / index ----------------------------------------------------

    def _u32(self) -> int:
        return struct.unpack(self._endian + "I", self._fh.read(4))[0]

    def _read_header(self) -> None:
        sig_bytes = self._fh.read(4)
        sig_le = struct.unpack("<I", sig_bytes)[0]
        sig_be = struct.unpack(">I", sig_bytes)[0]
        if sig_le == 0x1A412743:
            self._endian = "<"
        elif sig_be == 0x1A412743:
            self._endian = ">"
        else:
            raise ValueError(f"not a .2bit file: bad signature in {self.path}")
        version = self._u32()
        if version != 0:
            raise ValueError(f"unsupported .2bit version {version}")
        self._seq_count = self._u32()
        reserved = self._u32()
        if reserved != 0:
            raise ValueError(f"unexpected reserved field {reserved}")

    def _read_index(self) -> None:
        for _ in range(self._seq_count):
            name_size = self._fh.read(1)[0]
            name = self._fh.read(name_size).decode("ascii")
            offset = self._u32()
            self._index[name] = offset

    # --- public ------------------------------------------------------------

    def sequences(self) -> list[str]:
        return list(self._index)

    def has(self, name: str) -> bool:
        return name in self._index

    def fetch(self, name: str) -> tuple[str, int]:
        """Return (sequence_string, dna_size) for chromosome `name`.

        The returned string is uppercase ACGTN with lowercase letters for
        soft-masked (repeat) regions, matching UCSC's normal soft-masked FASTA
        convention.
        """
        if name not in self._index:
            raise KeyError(f"no such sequence: {name}")
        f = self._fh
        f.seek(self._index[name])
        dna_size = self._u32()

        n_block_count = self._u32()
        n_starts = list(struct.unpack(self._endian + f"{n_block_count}I",
                                      f.read(4 * n_block_count)))
        n_sizes  = list(struct.unpack(self._endian + f"{n_block_count}I",
                                      f.read(4 * n_block_count)))

        m_block_count = self._u32()
        m_starts = list(struct.unpack(self._endian + f"{m_block_count}I",
                                      f.read(4 * m_block_count)))
        m_sizes  = list(struct.unpack(self._endian + f"{m_block_count}I",
                                      f.read(4 * m_block_count)))

        _reserved = self._u32()

        packed_len = (dna_size + 3) // 4
        packed = f.read(packed_len)

        # Decode all bases as uppercase ACGT first; then overwrite N-blocks and
        # apply soft-mask (lowercase) over mask-blocks. Doing it as a bytearray
        # avoids quadratic string concat for large sequences.
        out = bytearray(dna_size)
        # Fast bulk decode via lookup table over single byte → 4 chars.
        # Building the lookup once amortizes the cost across all chromosomes.
        if not hasattr(TwoBitReader, "_byte_to_4bases"):
            tbl = []
            for b in range(256):
                tbl.append(bytes([
                    _CODES[(b >> 6) & 3],
                    _CODES[(b >> 4) & 3],
                    _CODES[(b >> 2) & 3],
                    _CODES[b & 3],
                ]))
            TwoBitReader._byte_to_4bases = tbl  # type: ignore[attr-defined]

        tbl = TwoBitReader._byte_to_4bases  # type: ignore[attr-defined]
        for i, b in enumerate(packed):
            out[i * 4 : i * 4 + 4] = tbl[b]
        del out[dna_size:]  # trim padding from the last byte

        # Apply N blocks (overwrite with 'N').
        for s, sz in zip(n_starts, n_sizes):
            out[s : s + sz] = b"N" * sz

        # Apply soft-mask: lowercase the masked ranges.
        # ord('a') - ord('A') == 32, so OR with 0x20 lowercases ASCII letters.
        for s, sz in zip(m_starts, m_sizes):
            for j in range(s, s + sz):
                out[j] |= 0x20

        return out.decode("ascii"), dna_size

    def write_fasta(self, name: str, out_path: Path, line_width: int = 50) -> int:
        """Write a single sequence as FASTA. Returns number of bases written."""
        seq, dna_size = self.fetch(name)
        with open(out_path, "w") as f:
            f.write(f">{name}\n")
            for i in range(0, dna_size, line_width):
                f.write(seq[i : i + line_width])
                f.write("\n")
        return dna_size

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> "TwoBitReader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def main(argv: list[str] | None = None) -> int:
    """CLI: list sequences in a 2bit, or extract one as FASTA.

    Examples:
        python -m twobit list /scratch2/.../hg38.2bit
        python -m twobit extract /scratch2/.../hg38.2bit chr19 hg38.chr19.fa
    """
    import argparse
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("list")
    pl.add_argument("twobit", type=Path)
    pe = sub.add_parser("extract")
    pe.add_argument("twobit", type=Path)
    pe.add_argument("chrom", type=str)
    pe.add_argument("out_fa", type=Path)
    args = p.parse_args(argv)

    with TwoBitReader(args.twobit) as t:
        if args.cmd == "list":
            for name in t.sequences():
                print(name)
            return 0
        if args.cmd == "extract":
            n = t.write_fasta(args.chrom, args.out_fa)
            print(f"{args.chrom}: {n:,} bp → {args.out_fa}")
            return 0
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

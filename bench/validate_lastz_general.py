#!/usr/bin/env python3
"""Self-consistency check for lastz --format=general output.

This is "check 2" from the T1 correctness audit (see README §11,
"T1: OpenMP batched parallel y-drop"). It verifies internal
consistency of each output row without re-running the aligner, and
flags anything that suggests the parallel path is producing
malformed output.

Per row we verify:
  - coordinates are in [0, seqlen] and start < end on both seqs
  - matches <= length (the "matches/length" column)
  - identity% is within 0.1 of 100 * matches / length
  - reported seq lengths match what the caller passes via name=len
  - score is strictly positive

We also count:
  - exact-row duplicates (real "did this come out twice?" cases)
  - rows sharing the same (name1, strand1, start1, end1,
    name2, strand2, start2, end2) tuple but with different
    score/match counts (the "two anchors trimmed to the same box"
    masking-off artifact; expected in variant 2a / parallel,
    expected 0 in baseline)

Usage:
    validate_lastz_general.py <output.general> [<name>=<len> ...]

Examples:
    validate_lastz_general.py /tmp/cc.par.out \
        chr10_40000000_50000000=10000000 \
        chr20_40000000_50000000=10000000

Exit code is 0 iff there are no per-row errors. Duplicates and
shared-bbox rows are reported but do not fail (they are properties
of masking-off, not parallelism bugs).
"""
import sys


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2

    path = sys.argv[1]
    seqlens = {}
    for kv in sys.argv[2:]:
        name, length = kv.split('=')
        seqlens[name] = int(length)

    n_total = 0
    n_errors = 0
    seen_exact = {}
    seen_bbox = {}

    with open(path) as f:
        # general format has a "#" header line
        header = f.readline()
        for line in f:
            n_total += 1
            line = line.rstrip('\n')
            cols = line.split('\t')
            try:
                (score, name1, strand1, len1, start1, end1,
                 name2, strand2, len2, start2, end2,
                 matches_over_len, identity_pct,
                 cov, cov_pct) = cols
            except ValueError:
                n_errors += 1
                if n_errors <= 5:
                    print(f"ERR line {n_total}: bad column count "
                          f"({len(cols)}, expected 15)")
                    print(f"      {line}")
                continue
            score = int(score)
            len1, len2 = int(len1), int(len2)
            start1, end1 = int(start1), int(end1)
            start2, end2 = int(start2), int(end2)
            m, L = matches_over_len.split('/')
            m, L = int(m), int(L)
            identity_pct = float(identity_pct.rstrip('%'))

            errs = []
            if start1 < 0 or end1 > len1 or start1 >= end1:
                errs.append(
                    f"bad seq1 coords [{start1},{end1}] vs len {len1}")
            if start2 < 0 or end2 > len2 or start2 >= end2:
                errs.append(
                    f"bad seq2 coords [{start2},{end2}] vs len {len2}")
            if name1 in seqlens and seqlens[name1] != len1:
                errs.append(
                    f"seq1 length mismatch: reported {len1}, "
                    f"expected {seqlens[name1]}")
            if name2 in seqlens and seqlens[name2] != len2:
                errs.append(
                    f"seq2 length mismatch: reported {len2}, "
                    f"expected {seqlens[name2]}")
            if m > L:
                errs.append(f"matches {m} > length {L}")
            if L <= 0:
                errs.append(f"non-positive length {L}")
            else:
                computed = 100.0 * m / L
                if abs(computed - identity_pct) > 0.1:
                    errs.append(
                        f"identity {identity_pct}% != "
                        f"computed {computed:.2f}%")
            if score <= 0:
                errs.append(f"non-positive score {score}")

            exact_key = line
            seen_exact[exact_key] = seen_exact.get(exact_key, 0) + 1

            bbox_key = (name1, strand1, start1, end1,
                        name2, strand2, start2, end2)
            seen_bbox[bbox_key] = seen_bbox.get(bbox_key, 0) + 1

            if errs:
                n_errors += 1
                if n_errors <= 5:
                    print(f"ERR line {n_total}: {'; '.join(errs)}")
                    print(f"      {line}")

    n_exact_dup_rows = sum(c for c in seen_exact.values() if c > 1)
    n_exact_dup_keys = sum(1 for c in seen_exact.values() if c > 1)
    n_bbox_collisions = sum(c for c in seen_bbox.values() if c > 1)
    n_bbox_collision_keys = sum(1 for c in seen_bbox.values() if c > 1)

    print(f"total alignments:                  {n_total}")
    print(f"per-row arithmetic errors:         {n_errors}")
    print(f"exact-row duplicates:              "
          f"{n_exact_dup_rows} rows across "
          f"{n_exact_dup_keys} unique duplicated rows")
    print(f"shared-bbox rows (masking-off art): "
          f"{n_bbox_collisions} rows across "
          f"{n_bbox_collision_keys} shared bboxes")
    return 0 if n_errors == 0 else 1


if __name__ == '__main__':
    sys.exit(main())

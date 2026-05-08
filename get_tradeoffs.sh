#!/usr/bin/env bash
# Reproduces the two tables in README §11 "Within-species regime change":
#
#   1. Headline numbers — hg38 vs T2T-CHM13 chr1:50M-60M (3 variants)
#   2. Bottleneck flip across regimes (bird-Z cross-species + 2 hg-vs-CHM13 rows)
#
# Runs the three hg-vs-CHM13 lastz_TS invocations from sweep_human1.sh, plus a
# bird-Z 10 Mbp baseline (cached if already present). Total wall: ~95 s if
# bird-Z is cached, ~155 s on first run. Single-thread, pinned to CPU 0.
#
# Prereqs (one-time):
#   make -C lastz build_lastz_substages
#   python3 bench/fetch_genomes.py --only hs1
#   # Plus chr1 extraction + 10 Mbp slicing for both hg-vs-CHM13 and bird-Z;
#   # see README §11 "Reproduce" blocks for the exact commands.

set -euo pipefail

# Refuse to be sourced — sourcing leaks `set -e` into your shell.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then : ; else
    echo "get_tradeoffs.sh: don't source me, run me as './get_tradeoffs.sh'" >&2
    return 1 2>/dev/null || exit 1
fi

# Always operate from the repo root, regardless of where we're invoked.
cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"

# ----------------------------------------------------------------------------
# Prereq check
# ----------------------------------------------------------------------------

LASTZ=$(realpath lastz/src/lastz_TS)
HG=$(realpath    bench/data_genomes/hg38.chr1_50_60mb.fa)
CHM13=$(realpath bench/data_genomes/hs1.chr1_50_60mb.fa)
BIRD_T=$(realpath bench/data_genomes/galGal6.chrZ_0_10mb.fa)
BIRD_Q=$(realpath bench/data_genomes/taeGut2.chrZ_0_10mb.fa)

[[ -x "$LASTZ" ]]   || { echo "missing $LASTZ — run 'make -C lastz build_lastz_substages'" >&2; exit 1; }
for f in "$HG" "$CHM13" "$BIRD_T" "$BIRD_Q"; do
    [[ -f "$f" ]] || { echo "missing $f — see README §11 'Reproduce' for slice commands" >&2; exit 1; }
done

# ----------------------------------------------------------------------------
# Run-one helper (mirrors sweep_human1.sh's run_one but parameterizes outdir)
# ----------------------------------------------------------------------------

run_one() {
    local outdir="$1" tag="$2"; shift 2
    local stem="${outdir}/${tag}"

    local rc=0
    set +e
    /usr/bin/time -v -o "${stem}.time.txt" \
        env LASTZ_STAGE_REPORT="${stem}.stage.txt" \
            taskset -c 0 \
            "$LASTZ" "$@" --format=maf- \
            > "${stem}.maf" 2>"${stem}.err"
    rc=$?
    set -e

    if [[ $rc -ne 0 ]]; then
        printf '  *** lastz exited %d for tag=%s ***\n' "$rc" "$tag" >&2
        printf '  --- %s.err (last 5 lines) ---\n' "$stem" >&2
        tail -n 5 "${stem}.err" | sed 's/^/    /' >&2
        return $rc
    fi
}

# ----------------------------------------------------------------------------
# 1) hg-vs-CHM13 sweep — three variants, always re-run
# ----------------------------------------------------------------------------

HG_OUT=bench/results/seed-sweep-hg-vs-chm13
mkdir -p "$HG_OUT"

echo "[1/4] hg-vs-CHM13 cross_species_default  (~47 s)..." >&2
run_one "$HG_OUT" cross_species_default "$HG" "$CHM13"

echo "[2/4] hg-vs-CHM13 bigger_seed_match15    (~27 s)..." >&2
run_one "$HG_OUT" bigger_seed_match15   "$HG" "$CHM13" --seed=match15

echo "[3/4] hg-vs-CHM13 hg_recommended         (~22 s)..." >&2
run_one "$HG_OUT" hg_recommended        "$HG" "$CHM13" --seed=match15 --step=20 --notransition

# ----------------------------------------------------------------------------
# 2) bird-Z 10 Mbp baseline — cached if already present
# ----------------------------------------------------------------------------

BIRD_OUT=bench/results/seed-sweep-bird10mb
mkdir -p "$BIRD_OUT"
BIRD_STAGE="${BIRD_OUT}/default_12of19.stage.txt"
BIRD_MAF="${BIRD_OUT}/default_12of19.maf"

if [[ -s "$BIRD_STAGE" ]] && [[ -s "$BIRD_MAF" ]]; then
    echo "[4/4] bird-Z default_12of19: cached, reusing." >&2
else
    echo "[4/4] bird-Z default_12of19: running (~62 s)..." >&2
    run_one "$BIRD_OUT" default_12of19 "$BIRD_T" "$BIRD_Q"
fi

# ----------------------------------------------------------------------------
# Field extractors
# ----------------------------------------------------------------------------

# Pull the trailing number from a "<key>: <value>" line in a stage report.
# Example: field foo.stage.txt "seed hit search:" → "18.235"
field() { grep "^$2" "$1" | awk '{print $NF}'; }

hsps() { grep -c '^a score=' "$1" || true; }

# Sum the alignment-size field (col 4 of the s-line right after each `a` line).
bp_total() {
    awk '/^a score=/ { getline; split($0, a, " "); sum += a[4] } END { print sum+0 }' "$1"
}

# Format <bp> as "X.XX Mbp" for table prettiness.
mbp() { awk -v n="$1" 'BEGIN { printf "%5.2f Mbp", n/1e6 }'; }

# Compute integer percentage a/b (rounded), as a string with no % sign.
pct() { awk -v a="$1" -v b="$2" 'BEGIN { if (b==0) print "?"; else printf "%.0f", 100*a/b }'; }

# ----------------------------------------------------------------------------
# Table 1 — Headline numbers
# ----------------------------------------------------------------------------

declare -a HEADLINE=(
    # label                      stage_file                                              maf_file
    "12of19          ${HG_OUT}/cross_species_default.stage.txt    ${HG_OUT}/cross_species_default.maf"
    "match15         ${HG_OUT}/bigger_seed_match15.stage.txt      ${HG_OUT}/bigger_seed_match15.maf"
    "match15+step20  ${HG_OUT}/hg_recommended.stage.txt           ${HG_OUT}/hg_recommended.maf"
)

echo
echo "Headline numbers — hg38 vs T2T-CHM13 chr1:50M-60M"
echo "================================================="
printf "%-16s  %8s  %10s  %12s  %6s  %12s\n" \
       "Variant" "Wall" "seed_hit" "gapped_ext" "HSPs" "bp_aligned"
echo "------------------------------------------------------------------------------"
for row in "${HEADLINE[@]}"; do
    read -r label stage maf <<< "$row"
    wall=$(field    "$stage" "total run time:")
    seedhit=$(field "$stage" "seed hit search:")
    gapped=$(field  "$stage" "gapped extension:")
    n=$(hsps "$maf")
    bp=$(bp_total "$maf")
    printf "%-16s  %7ss  %9ss  %11ss  %6d  %12s\n" \
           "$label" "$wall" "$seedhit" "$gapped" "$n" "$(mbp "$bp")"
done

# ----------------------------------------------------------------------------
# Table 2 — Bottleneck flip across regimes
# ----------------------------------------------------------------------------

echo
echo "Bottleneck flip across regimes (seed_hit vs. gapped_ext share of total)"
echo "======================================================================="
printf "%-50s  %-9s  %10s  %12s\n" \
       "Workload" "Identity" "seed_hit %" "gapped_ext %"
echo "------------------------------------------------------------------------------"

# Bird-Z cross-species default — the "seed_hit_search dominates" regime.
bird_total=$(field   "$BIRD_STAGE" "total run time:")
bird_seedhit=$(field "$BIRD_STAGE" "seed hit search:")
bird_gapped=$(field  "$BIRD_STAGE" "gapped extension:")
printf "%-50s  %-9s  %9s%%  %11s%%\n" \
       "Bird-Z 10 Mbp (cross-species default --seed=12of19)" \
       "~75%" \
       "$(pct "$bird_seedhit" "$bird_total")" \
       "$(pct "$bird_gapped"  "$bird_total")"

# hg-vs-CHM13 with the (wrong) cross-species defaults — already partial flip.
hg_def_stage="${HG_OUT}/cross_species_default.stage.txt"
hg_def_total=$(field   "$hg_def_stage" "total run time:")
hg_def_seedhit=$(field "$hg_def_stage" "seed hit search:")
hg_def_gapped=$(field  "$hg_def_stage" "gapped extension:")
printf "%-50s  %-9s  %9s%%  %11s%%\n" \
       "hg-vs-CHM13 (cross-species default — wrong choice)" \
       "~99.5%" \
       "$(pct "$hg_def_seedhit" "$hg_def_total")" \
       "$(pct "$hg_def_gapped"  "$hg_def_total")"

# hg-vs-CHM13 with match15+step20 — the "gapped_ext dominates" regime.
hg_rec_stage="${HG_OUT}/hg_recommended.stage.txt"
hg_rec_total=$(field   "$hg_rec_stage" "total run time:")
hg_rec_seedhit=$(field "$hg_rec_stage" "seed hit search:")
hg_rec_gapped=$(field  "$hg_rec_stage" "gapped extension:")
printf "%-50s  %-9s  %9s%%  %11s%%\n" \
       "hg-vs-CHM13 (match15 + step=20 + notransition)" \
       "~99.5%" \
       "$(pct "$hg_rec_seedhit" "$hg_rec_total")" \
       "$(pct "$hg_rec_gapped"  "$hg_rec_total")"

echo
echo "(raw artifacts in $HG_OUT and $BIRD_OUT)"

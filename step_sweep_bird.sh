#!/usr/bin/env bash
# Mirror of step_sweep.sh, but on the cross-species (~75% identity) workload:
# galGal6 vs taeGut2 chrZ:0-10Mbp.
#
# Uses --seed=12of19 (the lastz cross-species default) — match15 at 75 %
# identity would yield essentially no HSPs and turn the sweep into noise.
# step=1 is the cached baseline from get_tradeoffs.sh / earlier eu-stack work
# (bench/results/seed-sweep-bird10mb/default_12of19.*); we reuse it as the
# leftmost Pareto point and add step ∈ {20, 30, 50, 100}.
#
# Pinned to CPU 0. Total wall: ~3 min for the four new runs.

set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then : ; else
    echo "step_sweep_bird.sh: don't source me, run me as './step_sweep_bird.sh'" >&2
    return 1 2>/dev/null || exit 1
fi

cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"

LASTZ=$(realpath lastz/src/lastz_TS)
BIRD_T=$(realpath bench/data_genomes/galGal6.chrZ_0_10mb.fa)
BIRD_Q=$(realpath bench/data_genomes/taeGut2.chrZ_0_10mb.fa)

[[ -x "$LASTZ" ]]  || { echo "missing $LASTZ — run 'make -C lastz build_lastz_substages'" >&2; exit 1; }
[[ -f "$BIRD_T" ]] || { echo "missing $BIRD_T" >&2; exit 1; }
[[ -f "$BIRD_Q" ]] || { echo "missing $BIRD_Q" >&2; exit 1; }

OUTDIR=bench/results/step-sweep-bird10mb
mkdir -p "$OUTDIR"

# Cached step=1 anchor — produced earlier by the eu-stack / bird-Z baseline run.
BASELINE_STAGE=bench/results/seed-sweep-bird10mb/default_12of19.stage.txt
BASELINE_MAF=bench/results/seed-sweep-bird10mb/default_12of19.maf
[[ -s "$BASELINE_STAGE" && -s "$BASELINE_MAF" ]] || {
    echo "missing cached step=1 baseline. Run get_tradeoffs.sh first." >&2; exit 1; }

run_one() {
    local tag="$1"; shift
    local stem="${OUTDIR}/${tag}"
    local rc=0
    set +e
    /usr/bin/time -v -o "${stem}.time.txt" \
        env LASTZ_STAGE_REPORT="${stem}.stage.txt" \
            taskset -c 0 \
            "$LASTZ" "$BIRD_T" "$BIRD_Q" "$@" --format=maf- \
            > "${stem}.maf" 2>"${stem}.err"
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
        printf '  *** lastz exited %d for tag=%s ***\n' "$rc" "$tag" >&2
        tail -n 5 "${stem}.err" | sed 's/^/    /' >&2
        return $rc
    fi
}

STEPS=(20 30 50 100)
for s in "${STEPS[@]}"; do
    echo "[step=$s] running... (~30-60 s)" >&2
    run_one "step${s}" --seed=12of19 --step=$s
done

# ----------------------------------------------------------------------------
# Pretty-print
# ----------------------------------------------------------------------------

field() { grep "^$2" "$1" | awk '{print $NF}'; }
hsps()  { grep -c '^a score=' "$1" || true; }
bp_total() {
    awk '/^a score=/ { getline; split($0, a, " "); sum += a[4] } END { print sum+0 }' "$1"
}
mbp() { awk -v n="$1" 'BEGIN { printf "%5.2f Mbp", n/1e6 }'; }
pct() { awk -v a="$1" -v b="$2" 'BEGIN { if (b==0) print "?"; else printf "%.0f", 100*a/b }'; }

print_row() {
    local label="$1" stage="$2" maf="$3"
    local wall=$(field    "$stage" "total run time:")
    local seedhit=$(field "$stage" "seed hit search:")
    local gapped=$(field  "$stage" "gapped extension:")
    local n=$(hsps "$maf")
    local bp=$(bp_total "$maf")
    local seed_share=$(pct "$seedhit" "$wall")
    local gap_share=$(pct  "$gapped"  "$wall")
    printf "%-6s  %7ss  %9ss (%2s%%)  %9ss (%2s%%)  %6d  %12s\n" \
           "$label" "$wall" "$seedhit" "$seed_share" "$gapped" "$gap_share" "$n" "$(mbp "$bp")"
}

echo
echo "Pareto sweep — bird-Z 10 Mbp (galGal6 vs taeGut2 chrZ, ~75% identity)"
echo "                --seed=12of19 (cross-species default), varying --step"
echo "===================================================================================="
printf "%-6s  %8s  %14s        %14s       %6s  %12s\n" \
       "step" "Wall" "seed_hit (share)" "gapped_ext (share)" "HSPs" "bp_aligned"
echo "------------------------------------------------------------------------------------"

print_row "1"   "$BASELINE_STAGE"          "$BASELINE_MAF"
for s in "${STEPS[@]}"; do
    print_row "$s" "${OUTDIR}/step${s}.stage.txt" "${OUTDIR}/step${s}.maf"
done

echo
echo "(raw artifacts in $OUTDIR; step=1 reused from bench/results/seed-sweep-bird10mb)"

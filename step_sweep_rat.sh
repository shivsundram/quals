#!/usr/bin/env bash
# Mouse-rat sweep at the *threshold*: mm10.chr10[40-50M] vs rn6.chr20[40-50M].
# ~85 % identity in syntenic regions, dense synteny chosen by probing 5
# (mm.chr10 x rn.chr20) 10 Mbp pairs at --step=100 (see
# bench/results/step-sweep-rat10mb/probe/). The selected pair has the highest
# alignable density (~3.77 Mbp at step=100) — squarely in the gap-dominated
# zone but at ~85 % identity instead of rhesus's ~93 %. This pin-points where
# the seed -> gapped bottleneck crossover sits along the density axis.
#
# Pinned to CPU 0. Total wall: ~5 min for the five runs.

set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then : ; else
    echo "step_sweep_rat.sh: don't source me, run me as './step_sweep_rat.sh'" >&2
    return 1 2>/dev/null || exit 1
fi

cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"

LASTZ=$(realpath lastz/src/lastz_TS)
MM=$(realpath bench/data_genomes/mm10.chr10_40_50mb.fa)
RN=$(realpath bench/data_genomes/rn6.chr20_40_50mb.fa)

[[ -x "$LASTZ" ]]  || { echo "missing $LASTZ — run 'make -C lastz build_lastz_substages'" >&2; exit 1; }
[[ -f "$MM" ]]     || { echo "missing $MM"    >&2; exit 1; }
[[ -f "$RN" ]]     || { echo "missing $RN (run 'python3 bench/fetch_genomes.py --only rn6' + slice [40,50])" >&2; exit 1; }

OUTDIR=bench/results/step-sweep-rat10mb
mkdir -p "$OUTDIR"

run_one() {
    local tag="$1"; shift
    local stem="${OUTDIR}/${tag}"
    local rc=0
    set +e
    /usr/bin/time -v -o "${stem}.time.txt" \
        env LASTZ_STAGE_REPORT="${stem}.stage.txt" \
            taskset -c 0 \
            "$LASTZ" "$MM" "$RN" "$@" --format=maf- \
            > "${stem}.maf" 2>"${stem}.err"
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
        printf '  *** lastz exited %d for tag=%s ***\n' "$rc" "$tag" >&2
        tail -n 5 "${stem}.err" | sed 's/^/    /' >&2
        return $rc
    fi
}

STEPS=(1 20 30 50 100)
for s in "${STEPS[@]}"; do
    echo "[step=$s] running..." >&2
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

echo
echo "Pareto sweep — mm10.chr10[40-50M] x rn6.chr20[40-50M] (~85% identity, dense)"
echo "                --seed=12of19, varying --step"
echo "===================================================================================="
printf "%-6s  %8s  %14s        %14s       %6s  %12s\n" \
       "step" "Wall" "seed_hit (share)" "gapped_ext (share)" "HSPs" "bp_aligned"
echo "------------------------------------------------------------------------------------"

for s in "${STEPS[@]}"; do
    stage="${OUTDIR}/step${s}.stage.txt"
    maf="${OUTDIR}/step${s}.maf"
    wall=$(field    "$stage" "total run time:")
    seedhit=$(field "$stage" "seed hit search:")
    gapped=$(field  "$stage" "gapped extension:")
    n=$(hsps "$maf")
    bp=$(bp_total "$maf")
    seed_share=$(pct "$seedhit" "$wall")
    gap_share=$(pct  "$gapped"  "$wall")
    printf "%-6s  %7ss  %9ss (%2s%%)  %9ss (%2s%%)  %6d  %12s\n" \
           "$s" "$wall" "$seedhit" "$seed_share" "$gapped" "$gap_share" "$n" "$(mbp "$bp")"
done

echo
echo "(raw artifacts in $OUTDIR)"

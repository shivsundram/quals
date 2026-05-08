#!/usr/bin/env bash
# Sub-sweep along the --step axis at --seed=match15 --notransition on
# hg38 vs T2T-CHM13 chr1:50M-60M. Charts the sensitivity-vs-speed Pareto curve.
#
#   step=20   ← already in get_tradeoffs.sh as "hg_recommended"
#   step=30
#   step=50
#   step=100
#
# Single-thread, pinned to CPU 0. Total wall: ~75 s for these four runs.
#
# Prereqs:
#   make -C lastz build_lastz_substages
#   bench/data_genomes/hg38.chr1_50_60mb.fa
#   bench/data_genomes/hs1.chr1_50_60mb.fa

set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then : ; else
    echo "step_sweep.sh: don't source me, run me as './step_sweep.sh'" >&2
    return 1 2>/dev/null || exit 1
fi

cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"

LASTZ=$(realpath lastz/src/lastz_TS)
HG=$(realpath    bench/data_genomes/hg38.chr1_50_60mb.fa)
CHM13=$(realpath bench/data_genomes/hs1.chr1_50_60mb.fa)

[[ -x "$LASTZ" ]]  || { echo "missing $LASTZ — run 'make -C lastz build_lastz_substages'" >&2; exit 1; }
[[ -f "$HG" ]]     || { echo "missing $HG"    >&2; exit 1; }
[[ -f "$CHM13" ]]  || { echo "missing $CHM13" >&2; exit 1; }

OUTDIR=bench/results/step-sweep-hg-vs-chm13
mkdir -p "$OUTDIR"

run_one() {
    local tag="$1"; shift
    local stem="${OUTDIR}/${tag}"
    local rc=0
    set +e
    /usr/bin/time -v -o "${stem}.time.txt" \
        env LASTZ_STAGE_REPORT="${stem}.stage.txt" \
            taskset -c 0 \
            "$LASTZ" "$HG" "$CHM13" "$@" --format=maf- \
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
    echo "[step=$s] running..." >&2
    run_one "step${s}" --seed=match15 --step=$s --notransition
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

echo
echo "Pareto sweep — --seed=match15 --notransition, varying --step"
echo "============================================================="
printf "%-6s  %8s  %10s  %12s  %6s  %12s  %10s\n" \
       "step" "Wall" "seed_hit" "gapped_ext" "HSPs" "bp_aligned" "bp/HSP"
echo "------------------------------------------------------------------------------"
for s in "${STEPS[@]}"; do
    stage="${OUTDIR}/step${s}.stage.txt"
    maf="${OUTDIR}/step${s}.maf"
    wall=$(field    "$stage" "total run time:")
    seedhit=$(field "$stage" "seed hit search:")
    gapped=$(field  "$stage" "gapped extension:")
    n=$(hsps "$maf")
    bp=$(bp_total "$maf")
    mean=$(awk -v n="$n" -v bp="$bp" 'BEGIN { if (n==0) print "?"; else printf "%.0f", bp/n }')
    printf "%-6s  %7ss  %9ss  %11ss  %6d  %12s  %10s\n" \
           "$s" "$wall" "$seedhit" "$gapped" "$n" "$(mbp "$bp")" "$mean"
done

echo
echo "(raw artifacts in $OUTDIR)"

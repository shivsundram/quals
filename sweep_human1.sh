#!/usr/bin/env bash
# Within-species seed-config sweep on hg38 vs T2T-CHM13 chr1:50M-60M.
#
# Three lastz_TS runs with different seed configurations, plus a
# sensitivity-check loop using total-bp-aligned (not HSP count) — see
# README §11 "Within-species regime change" for the analysis these
# numbers feed.
#
# Prereqs (one-time):
#   make -C lastz build_lastz_substages
#   python3 bench/fetch_genomes.py --only hs1
#   python3 -m bench.twobit extract \
#       /scratch2/shiv1/lastz-bench-data/2bit/hg38.2bit chr1 \
#       /scratch2/shiv1/lastz-bench-data/slices/hg38.chr1.fa
#   ln -sf /scratch2/shiv1/lastz-bench-data/slices/hg38.chr1.fa \
#       bench/data_genomes/hg38.chr1.fa
#   python3 bench/slice_genome.py bench/data_genomes/hg38.chr1.fa \
#       bench/data_genomes/hg38.chr1_50_60mb.fa --start 50000000 \
#       --length 10000000 --rename hg38_chr1_50M_60M
#   python3 bench/slice_genome.py bench/data_genomes/hs1.chr1.fa \
#       bench/data_genomes/hs1.chr1_50_60mb.fa  --start 50000000 \
#       --length 10000000 --rename hs1_chr1_50M_60M
#
# Total wall: ~95 s on a single core (head node OK at this size).

set -euo pipefail

# Always operate relative to the repo root, regardless of where we're invoked.
cd "$(dirname "$(readlink -f "$0")")"

OUT=bench/results/seed-sweep-hg-vs-chm13
mkdir -p "$OUT"
cd "$OUT"

TARGET=$(realpath ../../../bench/data_genomes/hg38.chr1_50_60mb.fa)
QUERY=$(realpath  ../../../bench/data_genomes/hs1.chr1_50_60mb.fa)
LASTZ=$(realpath  ../../../lastz/src/lastz_TS)

# Fail fast with a useful message if any prereq is missing — without -e and
# these checks the script would silently produce empty .maf files.
[[ -x "$LASTZ" ]] || { echo "missing $LASTZ — run 'make -C lastz build_lastz_substages'" >&2; exit 1; }
[[ -f "$TARGET" ]] || { echo "missing $TARGET — see header for slice commands" >&2; exit 1; }
[[ -f "$QUERY"  ]] || { echo "missing $QUERY  — see header for slice commands" >&2; exit 1; }

run_one() {
    local tag="$1"; shift
    printf '\n=== %s : %s ===\n' "$tag" "$*"

    # Run lastz under /usr/bin/time, capturing stderr to a file. We disable
    # `set -e` for the duration of this command so we can inspect the exit
    # code and print a useful diagnostic before bailing — without this, a
    # bad --seed= or missing input fails silently because the .err file is
    # the only place the error appears.
    local rc=0
    set +e
    /usr/bin/time -v -o "${tag}.time.txt" \
        env LASTZ_STAGE_REPORT="${tag}.stage.txt" \
            taskset -c 0 \
            "$LASTZ" "$TARGET" "$QUERY" "$@" --format=maf- \
            > "${tag}.maf" 2>"${tag}.err"
    rc=$?
    set -e

    if [[ $rc -ne 0 ]]; then
        printf '  *** lastz exited %d for tag=%s ***\n' "$rc" "$tag" >&2
        printf '  --- %s.err (last 5 lines) ---\n' "$tag" >&2
        tail -n 5 "${tag}.err" | sed 's/^/    /' >&2
        printf '  hint: a 2-3 KB .maf with no "^a score=" lines is the\n' >&2
        printf '        --help text dumped to stdout when lastz rejects an arg.\n' >&2
        printf '        Common causes: --seed=matchN with N>15 (bit weight > 31),\n' >&2
        printf '        invalid --step= value, or missing/unreadable input file.\n' >&2
        return $rc
    fi

    local wall total seedhit
    wall=$(grep 'Elapsed (wall'    "${tag}.time.txt"  | awk '{print $NF}')
    total=$(grep '^total run time:'   "${tag}.stage.txt" | awk '{print $NF}')
    seedhit=$(grep '^seed hit search:' "${tag}.stage.txt" | awk '{print $NF}')
    printf '  wall=%s  total=%ss  seed_hit=%ss\n' "$wall" "$total" "$seedhit"
}

# NOTE on seed-size ceiling: lastz hardcodes maxSeedBitWeight=31 in
#   lastz/src/seeds.h:81. For contiguous seeds (bit_weight = 2*N) the
#   practical max is match15 (bit weight 30). match16+ get rejected at
#   startup with "N is not a valid word length", and lastz prints its
#   --help to stdout — which our run_one will capture into the .maf file
#   and report as "0 HSPs" if you don't notice the diagnostic. To try
#   bigger seeds, either:
#     (a) push --step further while keeping --seed=match15, or
#     (b) patch maxSeedBitWeight (and possibly maxSeedLen=31) and rebuild.
run_one cross_species_default
run_one bigger_seed_match15  --seed=match15
run_one hg_recommended       --seed=match15 --step=20 --notransition

# HSP count is misleading at sparser seeding (fragments merge into longer
# HSPs); total_bp is the honest sensitivity metric. See README §11 finding #2.
printf '\n=== sensitivity check (total_bp, not HSP count) ===\n'
for f in *.maf; do
    awk -v tag="${f%.maf}" '
        /^a score=/ { getline; split($0, a, " "); sum += a[4]; n++ }
        END { printf "%-25s HSPs=%d total_bp=%d mean=%d\n",
                     tag, n, sum, n ? sum/n : 0 }' "$f"
done

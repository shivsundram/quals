#!/bin/bash
#SBATCH --job-name=lastz-bench
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/shiv1/quals/bench/sbatch/logs/%x.%j.out
#SBATCH --error=/home/shiv1/quals/bench/sbatch/logs/%x.%j.err
#
# Run the lastz benchmark harness on a compute node.
#
# Required env vars (set via --export when sbatch'ing, or default below):
#   SUITE       — name of suite in workloads.py to run     (required)
#   REPS        — repetitions per cell                     (default: 1)
#   RUN_NAME    — bench output dir tag                     (default: <suite>-<jobid>)
#   PIN_CPU     — physical cpu id within allocation to pin (default: 0)
#
# We pin to a single CPU for clean timing. SLURM gives us 4 cpus_per_task so
# the OS scheduler has headroom; the benchmark itself runs on PIN_CPU only.

set -euo pipefail

SUITE="${SUITE:?must set SUITE=...}"
REPS="${REPS:-1}"
RUN_NAME="${RUN_NAME:-${SUITE}-${SLURM_JOB_ID}}"

# Pin to the FIRST cpu SLURM has actually given us (respects cgroup), unless
# the user overrode PIN_CPU explicitly.
ALLOC_CPUS="$(taskset -pc $$ | sed 's/.*: //')"
DEFAULT_PIN="$(echo "$ALLOC_CPUS" | cut -d, -f1 | cut -d- -f1)"
PIN_CPU="${PIN_CPU:-$DEFAULT_PIN}"

cd /home/shiv1/quals/bench

echo "==== job $SLURM_JOB_ID on $(hostname) ===="
echo "SUITE=$SUITE REPS=$REPS RUN_NAME=$RUN_NAME PIN_CPU=$PIN_CPU"
echo "Allocated CPUs:  $ALLOC_CPUS"
echo "Available mem:   $(free -h | awk '/Mem:/ {print $7}')"
echo "lastz binary:    $(ls -l ../lastz/src/lastz_T)"
echo "============================================"

python3 bench.py \
    --suite "$SUITE" \
    --reps  "$REPS" \
    --pin-cpu "$PIN_CPU" \
    --run-name "$RUN_NAME"

echo "==== done; results in bench/results/$RUN_NAME ===="

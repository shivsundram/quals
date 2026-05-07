#!/bin/bash
#SBATCH --job-name=lastz-profile-eustack
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/home/shiv1/quals/bench/sbatch/logs/%x.%j.out
#SBATCH --error=/home/shiv1/quals/bench/sbatch/logs/%x.%j.err
#
# Run bench/profile_eustack.py on a compute node — eu-stack PC sampling
# of a single lastz run. Typical use: profile the seed_hit_search stage
# at a manageable input size (e.g. the bird Z 10 Mbp slice, ~60-90 s wall).
#
# Required env vars (via --export when sbatch'ing):
#   TARGET     — target FASTA path (relative to repo root, or absolute)
#   QUERY      — query  FASTA path
#   OUT_DIR    — bench/results/<dir-name> for the profile output
#
# Optional:
#   LASTZ      — lastz binary               (default: ../lastz/src/lastz_T)
#   RATE_HZ    — sample rate (Hz)           (default: 100)
#   EXTRA_ARG  — single extra lastz arg     (default: --format=maf-)
#   PIN_CPU    — physical cpu to pin to     (default: first allocated)
#
# Examples:
#   # 10 Mbp bird-Z cross slice
#   sbatch --export=ALL,\
#       TARGET=bench/data_genomes/galGal6.chrZ_0_10mb.fa,\
#       QUERY=bench/data_genomes/taeGut2.chrZ_0_10mb.fa,\
#       OUT_DIR=bench/results/eustack-bird10mb-001 \
#       bench/sbatch/profile_eustack.sbatch
#
#   # Custom rate + extra lastz arg
#   sbatch --export=ALL,TARGET=...,QUERY=...,OUT_DIR=...,RATE_HZ=200,\
#       EXTRA_ARG=--nogapped \
#       bench/sbatch/profile_eustack.sbatch
#
# Notes on sizing:
#   bird Z cross 10 Mbp:  ~60-90 s wall  -> ~6000-9000 samples @ 100 Hz
#   bird Z cross full:    ~73 min wall   -> ~440k samples @ 100 Hz
#                                          (probably too many; bump --time
#                                           and lower RATE_HZ or use
#                                           --max-samples)
#
# We pin to a single CPU for clean wall time. SLURM gives us 4 cpus_per_task
# so the OS scheduler + the eu-stack sampler have headroom; the lastz
# subprocess itself runs on PIN_CPU only via taskset.

set -euo pipefail

TARGET="${TARGET:?must set TARGET=path/to/target.fa}"
QUERY="${QUERY:?must set QUERY=path/to/query.fa}"
OUT_DIR="${OUT_DIR:?must set OUT_DIR=bench/results/<run-name>}"
LASTZ="${LASTZ:-../lastz/src/lastz_T}"
RATE_HZ="${RATE_HZ:-100}"
EXTRA_ARG="${EXTRA_ARG:---format=maf-}"

# Pin to the FIRST cpu SLURM has actually given us (cgroup-aware).
ALLOC_CPUS="$(taskset -pc $$ | sed 's/.*: //')"
DEFAULT_PIN="$(echo "$ALLOC_CPUS" | cut -d, -f1 | cut -d- -f1)"
PIN_CPU="${PIN_CPU:-$DEFAULT_PIN}"

cd /home/shiv1/quals

echo "==== job $SLURM_JOB_ID on $(hostname) ===="
echo "TARGET=$TARGET"
echo "QUERY=$QUERY"
echo "OUT_DIR=$OUT_DIR"
echo "LASTZ=$LASTZ"
echo "RATE_HZ=$RATE_HZ"
echo "EXTRA_ARG=$EXTRA_ARG"
echo "PIN_CPU=$PIN_CPU  (allocated: $ALLOC_CPUS)"
echo "Time limit:  $(squeue -h -j $SLURM_JOB_ID -O TimeLimit)"
echo "============================================"

# We pin via taskset wrapping the python parent — both the python sampler
# AND the lastz child end up on PIN_CPU (eu-stack itself is a separate
# subprocess of the python script and inherits that affinity, which is fine
# because eu-stack's wall-clock fingerprint is tiny per sample).
taskset -c "$PIN_CPU" \
    python3 bench/profile_eustack.py \
        --lastz   "$LASTZ" \
        --target  "$TARGET" \
        --query   "$QUERY" \
        --out-dir "$OUT_DIR" \
        --rate-hz "$RATE_HZ" \
        --extra-arg="$EXTRA_ARG"

echo "==== done; results in $OUT_DIR ===="
echo "Quick look:"
echo "  cat $OUT_DIR/manifest.json"
echo "  cat $OUT_DIR/top_funcs.txt"
echo "  cat $OUT_DIR/top_lines.txt"
echo "  head -20 $OUT_DIR/folded.txt"

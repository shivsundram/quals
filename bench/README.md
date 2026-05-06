# LASTZ Benchmark Harness

Lightweight harness for running LASTZ against a configurable set of
**workloads** under a configurable set of **variants** (lastz CLI flag
combinations), capturing wall/user/sys time, peak RSS, page faults, and
context switches per run, and emitting structured results.

## Layout

```
bench/
  bench.py             # CLI runner
  workloads.py         # Workload + Variant + Suite registry
  gen_data.py          # Generate synthetic FASTA pairs (random DNA + mutated copy)
  fetch_real_data.py   # Pull curated real biological sequences from NCBI E-utilities
  data/                # Generated synthetic inputs (auto-created, idempotent)
  data_real/           # Cached real sequences (NCBI fetches; idempotent)
  results/<run>/       # Per-run output: manifest.json, summary.csv, runs/*.json
```

## Quick start

```bash
# 1) Build lastz first (once). Build BOTH the un-instrumented and stage-timer
#    versions; the harness prefers the timed binary when it's available.
make -C lastz                       # → lastz/src/lastz, lastz_D
make -C lastz build_lastz_timed     # → lastz/src/lastz_T  (with -DdbgTiming
                                    #   and -DdbgTimingGappedExtend)

# 2) See what's defined.
python3 bench/bench.py --list

# 3) Run the smoke suite (5 cells x 3 reps, ~10s total on the test_data inputs).
python3 bench/bench.py --suite smoke --reps 3 --pin-cpu 0

# 4) Run the synthetic scaling sweep (10 kbp .. 10 Mbp, ~70s total).
python3 bench/bench.py --suite scaling --reps 3 --pin-cpu 0 \
    --run-name scaling-baseline

# 5) One-off: a single (workload, variant) pair with custom reps.
python3 bench/bench.py --workload synth.uniform_1000kb --variant nogapped \
    --reps 5 --pin-cpu 0

# 6) Force the un-instrumented binary (useful for measuring the overhead of
#    stage timing itself, or for clean wall-time numbers).
python3 bench/bench.py --suite scaling --no-stage-timers --reps 3

# 7) Fetch real biological sequences (β-actin mRNA from NCBI for human,
#    chimp, and fruit fly) and run the 3x3 cross-species matrix.
python3 bench/fetch_real_data.py
python3 bench/bench.py --suite actb_matrix --reps 3 --keep-output
```

## Concepts

- **Workload**: a target/query input pair plus any required lastz flags. Two
  kinds:
    - `static.*` — points at files in `lastz/test_data/`.
    - `synth.uniform_<N>kb` — generated on demand by `gen_data.py`. Outputs
      are deterministic in `(length, divergence, indel_rate, seed)`, so re-runs
      reuse cached files (verified via SHA256 in each run's JSON).
- **Variant**: extra lastz flags layered on top, used to isolate stage costs:
    - `default` — all lastz defaults (seed-and-extend + gapped extension).
    - `nogapped` / `ungapped_only` — skip gapped extension.
    - `chain` — enable chaining stage.
    - `step20` — sparser seeds (lower sensitivity, faster).
- **Suite**: a named list of `(workload, variant)` cells.
- **Rep**: each cell is run `--reps` times after `--warmup` discarded reps.
  Aggregation reports min/median/max wall time and median RSS / page faults /
  output bytes.

## Output

For each run (one invocation of `bench.py`):

```
results/<run-name>/
  manifest.json              # env (host, cpu count, lastz path+version,
                             #   lastz_instrumented), argv, cells, started_at
  runs/<wl>__<var>__<rep>.json
                             # per-rep details: full cmd, exit code, all timing
                             #   + memory metrics, lastz stderr tail, input
                             #   sha256s (workload_meta), per-stage seconds
                             #   (when run with the instrumented binary)
  runs/<wl>__<var>__<rep>.stage.txt
                             # raw stage report from lastz (instrumented runs
                             #   only) — written by lastz itself via the
                             #   LASTZ_STAGE_REPORT env var
  summary.csv                # one row per (workload, variant) cell
```

`summary.csv` core columns:

| column | meaning |
|---|---|
| `wall_s_min/median/max` | wall-clock seconds, aggregated across reps |
| `user_s_median`, `sys_s_median` | CPU time decomposition |
| `cpu_pct_median` | `time -v`'s "% of CPU" — useful sanity check for pinning |
| `max_rss_kb_median` | peak resident set size |
| `page_faults_major/minor_median` | I/O vs in-memory faults |
| `output_bytes_median` | size of lastz output (proxy for # alignments) |
| `ok` | 1 if all reps in cell exited cleanly, else 0 |

When run with the instrumented binary (`lastz_T`), additional per-stage
median columns are appended (only those that have data appear):

| column | meaning |
|---|---|
| `total_run_time_s_median` | lastz's own end-to-end wall (matches `wall_s_median` modulo `time -v` overhead) |
| `sequence_1_io_s_median` / `sequence_2_io_s_median` | reading target / query |
| `seed_position_table_s_median` | building the target k-mer index |
| `seed_hit_search_s_median` | finding seed hits + ungapped extension to HSPs |
| `chaining_s_median` | chaining HSPs (only with `--chain`) |
| `gapped_extension_s_median` | gapped extension (y-drop banded SW DP) |
| `interpolation_s_median` | tweener / inner alignment between gapped pieces |
| `output_s_median` | formatting + writing alignments |
| `total_query_time_s_median` | per-query work, summed across queries |
| `ge_ydrop_align_s_median` | top-level y-drop DP driver inside gapped extension |
| `ge_ydrop_one_sided_align_s_median` | the actual inner DP loop (the GPU target) |
| `ge_update_lr_bounds_s_median`, `ge_update_active_segs_s_median`, etc. | bookkeeping inside gapped extension |

## Reproducibility

- Synthetic inputs are seeded; SHA256s of target+query are stored in each
  per-rep JSON so we can confirm two runs saw identical bytes.
- `--pin-cpu K` pins the lastz process to one core via `taskset -c K`.
- `--warmup N` runs N rep(s) first whose results are discarded — fills the
  page cache so measured reps see steady-state I/O.
- The lastz binary path and `--version` string are recorded in `manifest.json`.

## Adding a workload

Edit `workloads.py`:

```python
_register(Workload(
    id="static.my_input",
    description="...",
    resolve=_static("path/in/lastz/test_data/target.fa",
                    "path/in/lastz/test_data/query.fa"),
    tags=("smoke",),
    target_actions=("multiple",),    # if target has >1 sequence
    lastz_args=("--format=maf",),    # workload-specific flags
))
```

Synthetic workloads are similar but call `_synthetic(name, target_len, ...)`.

## Adding a variant

```python
VARIANTS["my_variant"] = Variant(
    id="my_variant",
    description="...",
    extra_args=["--step=10", "--seed=match12"],
)
```

## Stage timing — how it works

The instrumented binary (`lastz_T`) is built with `-DdbgTiming` and
`-DdbgTimingGappedExtend`, which switch on accumulator-style timers that
are already present in `lastz.c` and `gapped_extend.c`. These wrap each
stage with `clock_gettime`-style start/stop pairs and accumulate elapsed
microseconds; on shutdown lastz prints a per-stage report.

Two small additions on top of the upstream code (in `lastz.c`):

1. The report is bracketed with `===STAGE_TIMING_BEGIN===` /
   `===STAGE_TIMING_END===` markers, so a parser can find it
   deterministically inside stderr.
2. If the env var `LASTZ_STAGE_REPORT=<path>` is set, the report is
   written to that file instead of stderr (so the harness doesn't have to
   filter through the rest of lastz's stderr).

Correctness: the timing instrumentation is purely bookkeeping (no control
flow changes), and was verified to produce byte-identical alignments to
the un-instrumented build on the bundled `make test` reference.

Overhead: the timers fire only at coarse stage boundaries (a few hundred
calls for typical workloads), so the cost is sub-millisecond and well
below run-to-run variance. If you want to confirm, run the same suite
with `--no-stage-timers` and compare wall medians.

## Known caveats

- For target sizes above ~1 Mbp, lastz emits warnings about truncating
  alignments due to traceback memory limits. The harness captures this in
  `stderr_tail` of each run's JSON; bump traceback with
  `--allocate:traceback=<bytes>` (via a variant) when running large sizes
  for clean numbers.
- `perf` is not usable on this kernel (paranoid=3, missing linux-tools for
  5.4.0-216), so we rely on `/usr/bin/time -v` plus the in-binary stage
  timers. Hardware counters (cache misses, IPC, branch-pred) will need
  either (a) sysctl change, (b) a different host, or (c) further
  source-level instrumentation in lastz itself.
- Single-threaded only — lastz is itself single-threaded, and the harness
  runs cells sequentially.

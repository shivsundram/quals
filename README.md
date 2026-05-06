# LASTZ → GPU acceleration (quals project)

End-to-end workspace for profiling [LASTZ](https://lastz.github.io/lastz/)
(a single-threaded CPU pairwise DNA aligner) and figuring out which stages
to GPU-accelerate. The longer-term goal is to reproduce / improve on
[SegAlign](https://github.com/gsneha26/SegAlign), the SC20 GPU LASTZ
accelerator, and to think about which serial stages might justify custom
hardware.

This README is the **single source of truth for how to recreate the
environment from scratch** on any machine — what to build, where data
comes from, how to run benchmarks, and what we already know.

---

## Table of contents

1. [System requirements](#1-system-requirements)
2. [Repo layout](#2-repo-layout)
3. [One-shot bootstrap](#3-one-shot-bootstrap-fresh-machine)
4. [Building lastz](#4-building-lastz)
5. [Data layer](#5-data-layer)
   - [5.1 Synthetic FASTAs (`gen_data.py`)](#51-synthetic-fastas-gen_datapy)
   - [5.2 Real mRNAs from NCBI (`fetch_real_data.py`)](#52-real-mrnas-from-ncbi-fetch_real_datapy)
   - [5.3 Reference genomes from UCSC (`fetch_genomes.py` + `twobit.py`)](#53-reference-genomes-from-ucsc-fetch_genomespy--twobitpy)
6. [Workload registry (`workloads.py`)](#6-workload-registry-workloadspy)
7. [Benchmark harness (`bench.py`)](#7-benchmark-harness-benchpy)
8. [Stage-timing instrumentation](#8-stage-timing-instrumentation)
9. [Running on SLURM](#9-running-on-slurm)
10. [What we know so far](#10-what-we-know-so-far)
11. [Open questions / next steps](#11-open-questions--next-steps)
12. [Gotchas, lessons learned, and "don't do that"](#12-gotchas-lessons-learned-and-dont-do-that)

---

## 1. System requirements

Tested on Ubuntu 18.04-era Linux (kernel 5.4, glibc 2.27) plus a SLURM
cluster with a `gpu` partition (Pascal cards, 4×/node). Anything modern
should work; the things that actually matter:

| Need | Why | How to verify |
|---|---|---|
| `gcc` ≥ 4.8 | builds lastz | `gcc --version` |
| `python3` ≥ 3.8 | bench harness, fetchers | `python3 --version` |
| `make` | builds lastz | `make --version` |
| `/usr/bin/time` (GNU, not the bash builtin) | bench harness reads `time -v` output | `/usr/bin/time --version` |
| `taskset` | CPU pinning | `taskset --version` |
| `curl` or `wget` | UCSC + NCBI downloads | already on most boxes |
| Network egress to `eutils.ncbi.nlm.nih.gov` and `hgdownload.soe.ucsc.edu` | fetching data | `curl -I https://hgdownload.soe.ucsc.edu` |
| Scratch directory at `/scratch2/shiv1/lastz-bench-data/` (or change the path) | genomes are several GB; don't put them under `$HOME` | `mkdir -p /scratch2/shiv1/lastz-bench-data` |
| SLURM with a `gpu` partition (optional, for the chromosome-scale runs) | head node is for editing, not running | `sinfo` |

No Python packages required beyond stdlib. No Conda env. No Docker. No
biopython.

---

## 2. Repo layout

```
quals/
├── README.md                          # this file
├── lastz/                             # upstream LASTZ source, slightly patched
│   ├── Makefile                       # +build_lastz_timed target
│   ├── src/
│   │   ├── lastz.c                    # +stage-timing markers + LASTZ_STAGE_REPORT env var
│   │   ├── Makefile                   # +%_T.o rule + lastz_T target (stage timers ON)
│   │   ├── lastz                      # built: production binary
│   │   ├── lastz_D                    # built: debug-symbols binary (upstream)
│   │   └── lastz_T                    # built: stage-timer-instrumented binary  ← we use this
│   └── test_data/                     # bundled smoke-test FASTAs (~kbp scale)
│
└── bench/
    ├── README.md                      # detail on the harness internals
    ├── bench.py                       # main runner: (workload × variant × rep) → metrics
    ├── workloads.py                   # workload + variant + suite registry
    ├── gen_data.py                    # synthetic FASTA generator
    ├── fetch_real_data.py             # NCBI mRNA fetcher
    ├── fetch_genomes.py               # UCSC 2bit fetcher + chromosome slicer
    ├── twobit.py                      # in-tree .2bit reader (replaces UCSC twoBitToFa)
    ├── data/                          # synthetic FASTAs (auto-created)
    ├── data_real/                     # NCBI fetches (auto-created)
    ├── data_genomes/                  # symlinks → /scratch2/.../2bit and slices
    ├── results/<run-name>/            # per-run output directory
    └── sbatch/
        ├── run_bench.sbatch           # SLURM wrapper — submits a suite to the gpu partition
        └── logs/                      # sbatch stdout/stderr per job
```

Plus, off the workspace:

```
/scratch2/shiv1/lastz-bench-data/      # large data lives here; symlinked into bench/data_genomes
├── 2bit/                              # hg38, mm10, galGal6, taeGut2 + md5 manifests
├── slices/                            # per-chromosome FASTAs extracted from each 2bit
└── tools/                             # downloaded UCSC binaries (e.g. twoBitToFa — not used, see §12)
```

---

## 3. One-shot bootstrap (fresh machine)

```bash
# (a) Get the workspace. lastz is a submodule pinned at our patched fork.
git clone --recurse-submodules git@github.com:shivsundram/quals.git
cd quals

# If you forgot --recurse-submodules:
#     git submodule update --init --recursive

# (b) Build both lastz binaries (production + stage-timer instrumented).
make -C lastz                           # → lastz/src/lastz, lastz_D
make -C lastz build_lastz_timed         # → lastz/src/lastz_T

# (c) Make scratch and create symlink target.
mkdir -p /scratch2/shiv1/lastz-bench-data/{2bit,slices,tools}

# (d) Pull the real data we benchmark against.
python3 bench/fetch_real_data.py        # ~few KB; mRNAs into bench/data_real/
python3 bench/fetch_genomes.py          # ~5 GB; genomes into scratch + symlinks

# (e) Smoke-test the harness on bundled tiny inputs.
python3 bench/bench.py --suite smoke --reps 3 --pin-cpu 0

# (f) Synthetic scaling sweep (10 kbp .. 10 Mbp; ~70 s on a fast core).
python3 bench/bench.py --suite scaling --reps 3 --pin-cpu 0 \
    --run-name scaling-bootstrap

# (g) The big real-data run goes to the GPU partition via SLURM:
sbatch --export=ALL,SUITE=bird_z_cross,REPS=1,RUN_NAME=bird_z_cross-$(date +%Y-%m-%d) \
       bench/sbatch/run_bench.sbatch
squeue -u $USER                          # watch it
tail -f bench/sbatch/logs/lastz-bench.<jobid>.err
```

After that, results are under `bench/results/<run-name>/summary.csv`.

---

## 4. Building lastz

LASTZ ships a plain Makefile and uses no external libraries. We have
**three** binaries we care about:

| Binary | How | Purpose |
|---|---|---|
| `lastz/src/lastz` | `make -C lastz` | Production, no debug overhead. The "fair" target for any timing. |
| `lastz/src/lastz_D` | (built by upstream Makefile alongside `lastz`) | Same code, with `-g` debug symbols. |
| `lastz/src/lastz_T` | `make -C lastz build_lastz_timed` | **Stage-timer instrumented** (`-DdbgTiming -DdbgTimingGappedExtend`). The harness uses this by default — see §8. |

**What we patched in the upstream tree** (very small surface area):

- `lastz/src/Makefile` — added a `stageTimers=ON` knob, a parallel
  `%_T.o` object rule that always defines the timing macros, a
  `lastz_T` link target, and a matching `clean_builds` line.
- `lastz/Makefile` (top-level) — added a `build_lastz_timed` convenience
  target that just calls `make -C src lastz_T`.
- `lastz/src/lastz.c` — three changes, all behind `#ifdef dbgTiming`:
  1. The `dbg_timing_report` macro now takes a `FILE*` argument.
  2. The end-of-run report is bracketed with
     `===STAGE_TIMING_BEGIN===` / `===STAGE_TIMING_END===` markers so
     the harness can find it deterministically.
  3. If the env var `LASTZ_STAGE_REPORT=<path>` is set, the report is
     written to that file instead of stderr.

The instrumentation only fires at coarse stage boundaries (a few hundred
calls per typical run), so the runtime overhead is sub-millisecond and
well below run-to-run variance — verified by running both binaries on
the bundled smoke test and getting byte-identical alignments.

To rebuild from scratch:

```bash
make -C lastz/src clean_builds
make -C lastz                           # production
make -C lastz build_lastz_timed         # instrumented
```

To verify the binaries:

```bash
lastz/src/lastz   --version | head -1
lastz/src/lastz_T --version | head -1   # should print same version
make -C lastz test                      # upstream smoke test
```

---

## 5. Data layer

We have three independent data sources, used for three different things:

| Source | Script | Output | Size | What it's for |
|---|---|---|---|---|
| **Synthetic** | `gen_data.py` | `bench/data/uniform_<N>kb.{target,query}.fa` | KB → tens of MB | Scaling sweeps with controlled divergence |
| **Real mRNA from NCBI** | `fetch_real_data.py` | `bench/data_real/actb_*.fa` | KB | Cross-species sensitivity matrix; tiny, fast iteration |
| **Reference genomes from UCSC** | `fetch_genomes.py` + `twobit.py` | `/scratch2/.../2bit/*.2bit` + `slices/*.fa`, symlinked into `bench/data_genomes/` | GB | Real chromosome-scale benchmarks (SegAlign-comparable) |

### 5.1 Synthetic FASTAs (`gen_data.py`)

Generate a random target, then derive a query by applying controlled
substitutions and short indels. Used for scaling experiments where we
want to vary input size while holding everything else constant.

```bash
# Direct CLI (for one-off generation)
python3 bench/gen_data.py --name uniform_5mbp --length 5_000_000 \
    --sub-rate 0.05 --indel-rate 0.01 --seed 0xC0FFEE \
    --out-dir bench/data/

# But normally you don't call this — workloads.py does it lazily.
# Just run e.g. `bench.py --workload synth.uniform_1000kb` and the
# files materialize in bench/data/ on first use.
```

**Determinism is the point.** Output bytes depend only on
`(name, length, sub_rate, indel_rate, seed)`. `bench.py` records SHA256
of every input it fed to lastz, so two runs are byte-identical or
they're not.

Default sweep sizes (in `workloads.py`):
`synth.uniform_{10,100,1000,10000}kb` — span 10 kbp → 10 Mbp.

### 5.2 Real mRNAs from NCBI (`fetch_real_data.py`)

Pulls a small panel of homologous gene sequences via NCBI E-utilities
(no biopython dependency, just stdlib `urllib`). Cached on disk;
re-runs without `--force` are no-ops.

Currently fetched: **β-actin (ACTB) mRNA across human / chimp / fly**.

| Species | Accession | Length | Why |
|---|---|---|---|
| Human (Homo sapiens) | `NM_001101.5` | ~1.8 kb | Reference |
| Chimp (Pan troglodytes) | `NM_001009945.2` | ~1.8 kb | 1:1 ortholog of human ACTB |
| Fly (Drosophila melanogaster, Act5C) | `NM_001014725.2` | ~1.9 kb | Canonical cytoplasmic actin in fly; deeply diverged from vertebrate ACTB |

These are comically small (kb), but useful for two things:
- Fast iteration on output-format / parsing changes (whole 3×3 matrix
  runs in <1 s).
- Showing how lastz's defaults handle different evolutionary distances:
  human↔chimp gives one full-length alignment, human/chimp↔fly gives
  ~zero hits at default sensitivity (tuned for >70% identity).

```bash
python3 bench/fetch_real_data.py            # idempotent
python3 bench/fetch_real_data.py --force    # re-download
ls bench/data_real/                         # actb_human.fa, actb_chimp.fa, actb_fly.fa
```

If you ever extend this to other genes, add a new `Target` to the
`PANELS["actb"]` list (or make a new panel) — see the file's
`accessions=[primary, fallback, ...]` pattern, which protects you from
RefSeq accession bumps.

### 5.3 Reference genomes from UCSC (`fetch_genomes.py` + `twobit.py`)

Downloads the four assemblies SegAlign benchmarks against, into scratch:

| Assembly | Species | Size | Slice we extract |
|---|---|---|---|
| `hg38` | Human (Homo sapiens) | 3.08 Gbp | `chr19` (~58 Mbp) |
| `mm10` | Mouse (Mus musculus) | 2.72 Gbp | `chr10` (~130 Mbp) |
| `galGal6` | Chicken (Gallus gallus) | 1.05 Gbp | `chrZ` (~82 Mbp) |
| `taeGut2` | Zebra Finch (Taeniopygia guttata) | 1.02 Gbp | `chrZ` (~73 Mbp) |

Why these four: they're exactly what SegAlign uses in their paper, so
our wall-times are directly comparable.

The slices were chosen for evolutionary cross-pairs:

- `hg38.chr19` ↔ `mm10.chr10` — mostly syntenic mammal pair (~90 Mya
  divergence), heavy repeat content.
- `galGal6.chrZ` ↔ `taeGut2.chrZ` — sex-chromosome bird pair (~100 Mya
  divergence), much lower repeat fraction. **First real chromosome-scale
  benchmark we run.**

```bash
python3 bench/fetch_genomes.py                  # download + verify md5 + slice
python3 bench/fetch_genomes.py --force          # re-download everything
python3 bench/fetch_genomes.py --skip-slices    # 2bit only
```

After this:

```
/scratch2/shiv1/lastz-bench-data/
├── 2bit/
│   ├── hg38.2bit  (~830 MB)
│   ├── mm10.2bit  (~720 MB)
│   ├── galGal6.2bit  (~280 MB)
│   ├── taeGut2.2bit  (~270 MB)
│   └── md5sum.<asm>.txt   # verified against UCSC's published manifests
└── slices/
    ├── hg38.chr19.fa     (58 MB)
    ├── mm10.chr10.fa    (130 MB)
    ├── galGal6.chrZ.fa   (84 MB)
    └── taeGut2.chrZ.fa   (74 MB)

bench/data_genomes/
├── hg38.2bit          → /scratch2/.../2bit/hg38.2bit
├── hg38.chr19.fa      → /scratch2/.../slices/hg38.chr19.fa
├── ... etc ...
```

#### Why the in-tree `twobit.py` reader

UCSC ships a `twoBitToFa` binary for extracting per-chromosome FASTAs
from .2bit files. **It doesn't run on this host** — it requires
`GLIBC_2.33+` which our 5.4 / Ubuntu-18.04-era kernel doesn't have.
Rather than ship a libc, we wrote a minimal pure-Python 2bit reader
(`bench/twobit.py`, ~200 lines) that:

- Parses the .2bit header + per-sequence index.
- Decodes 2-bit packed DNA (T=00, C=01, A=10, G=11).
- Honors **N-blocks** (assembly gaps → `N`).
- Honors **mask-blocks** (soft-masked repeat ranges → lowercase).
- Writes 60-col-wrapped FASTA matching what UCSC's twoBitToFa would
  produce.

It's not optimised — extracting chr10 of mm10 takes ~30 s — but it
runs on any Python 3 and we only do it once per assembly.

You can also use `twobit.py` as a CLI:

```bash
python3 bench/twobit.py /scratch2/.../2bit/hg38.2bit list
python3 bench/twobit.py /scratch2/.../2bit/hg38.2bit fetch chr22 > chr22.fa
```

---

## 6. Workload registry (`workloads.py`)

Single Python file that defines:

- **Workload** — a target/query input pair, plus any required lastz
  flags (e.g. `target_actions=("multiple",)` for multi-sequence FASTAs,
  `lastz_args=("--format=maf",)` to override default lav format).
- **Variant** — extra lastz CLI flags layered on top, used to isolate
  stage costs.
- **Suite** — a named list of `(workload, variant)` cells.

### Current workloads

```
static.pseudocat_vs_pseudopig          tiny smoke test (kbp)
static.pseudopig_self                  multi-seq smoke test
static.reads_vs_pseudopig              short-reads-vs-genome shape

real.actb_<sp1>_vs_<sp2>               9-cell mRNA matrix (~1.8 kbp each)
                                       sp ∈ {human, chimp, fly}

real.galGal6_chrZ_self                 ~82 Mbp self-vs-self
real.galGal6_vs_taeGut2_chrZ           ~82 × ~73 Mbp bird-vs-bird   ← headline real workload
real.hg38_vs_mm10                      ~58 × ~130 Mbp mammal pair (heavier; longer)

synth.uniform_{10,100,1000,10000}kb    scaling sweep, 10 kbp → 10 Mbp
```

### Current variants

```
default            lastz with all defaults (seed-and-extend + gapped, no chaining)
segalign_default   default + --format=maf-   ← byte-comparable to SegAlign output
nogapped           --nogapped — isolates seed + ungapped extension
chain              --chain   — adds chaining stage on top of default
step20             --step=20 — sparser seeds; lower sensitivity, faster
```

**Important context on `segalign_default`:** SegAlign uses lastz's
defaults across the board (verified against their `src/main.cpp`
[gsneha26/SegAlign main, lines 75–98](https://github.com/gsneha26/SegAlign/blob/main/src/main.cpp)).
The only practical difference is output format (`maf-` instead of the
default `lav`). So our `segalign_default` variant is the canonical
"comparable to SegAlign" run.

### Current suites

```
smoke                       quick sanity sweep, ~10 s
scaling                     synthetic scaling, default flags only
scaling_nogapped            same but --nogapped (subtract for gapped cost)
actb_matrix                 9-cell ACTB cross-species matrix at default sensitivity
bird_z_self                 galGal6.chrZ self-vs-self, segalign_default
bird_z_self_breakdown       same workload, default + nogapped
bird_z_cross                galGal6.chrZ vs taeGut2.chrZ, segalign_default  ← key real run
bird_z_cross_breakdown      same workload, default + nogapped
```

### Adding a workload

```python
_register(Workload(
    id="real.my_input",
    description="...",
    resolve=_real(target_path=DATA_GENOMES / "x.fa",
                  query_path=DATA_GENOMES / "y.fa"),
    tags=("real", "genome"),
    target_actions=("multiple",),    # if target FASTA has >1 sequence
    lastz_args=("--format=maf",),    # if you need a workload-wide flag
))
```

### Adding a variant

```python
VARIANTS["my_variant"] = Variant(
    id="my_variant",
    description="...",
    extra_args=["--ydrop=15000", "--gappedthresh=2200"],
)
```

### Adding a suite

```python
SUITES["my_suite"] = [
    ("real.galGal6_vs_taeGut2_chrZ", "default"),
    ("real.galGal6_vs_taeGut2_chrZ", "nogapped"),
]
```

---

## 7. Benchmark harness (`bench.py`)

Runs `(workload × variant × rep)` cells under `/usr/bin/time -v`,
captures wall / user / sys / RSS / page-faults / context-switches, parses
the per-stage timing report from `lastz_T`, and writes structured
output to `bench/results/<run-name>/`.

### CLI

```bash
# List everything defined.
python3 bench/bench.py --list

# Run a suite.
python3 bench/bench.py --suite scaling --reps 3 --pin-cpu 0 \
    --run-name scaling-baseline-2026-05-06

# Run one cell.
python3 bench/bench.py --workload synth.uniform_1000kb --variant nogapped \
    --reps 5 --pin-cpu 0

# Force the un-instrumented binary (measure stage-timer overhead).
python3 bench/bench.py --suite scaling --no-stage-timers --reps 3

# Use a different lastz binary explicitly.
python3 bench/bench.py --lastz /path/to/other/lastz --workload ...
```

Useful flags:

| Flag | What |
|---|---|
| `--reps N` | repetitions per cell (median is reported in summary) |
| `--warmup N` | run N rep(s) first whose results are discarded — fills page cache |
| `--pin-cpu K` | pin lastz to CPU K via `taskset` (cgroup-aware on SLURM) |
| `--keep-output` | keep the lastz alignment output file (default: discard after measuring its size) |
| `--no-stage-timers` | force the un-instrumented `lastz` binary |
| `--lastz <path>` | override which lastz binary to use |
| `--run-name NAME` | name the output directory; defaults to a timestamp |

### Output layout

```
bench/results/<run-name>/
├── manifest.json            # env (host, cpus, lastz path+version, instrumented y/n),
│                            #   argv, cells, timestamps
├── runs/
│   ├── <wl>__<var>__<rep>.json   # per-rep details: full cmd, exit code,
│   │                              #   all timing+memory metrics, lastz stderr tail,
│   │                              #   workload meta (input SHA256s),
│   │                              #   per-stage seconds (instrumented runs)
│   ├── <wl>__<var>__<rep>.stage.txt   # raw stage report from lastz_T
│   └── <wl>__<var>__<rep>.out         # alignment output (if --keep-output)
└── summary.csv              # one row per (workload, variant) cell, aggregated
```

`summary.csv` columns include:

- `wall_s_{min,median,max}`, `user_s_median`, `sys_s_median`, `cpu_pct_median`
- `max_rss_kb_median`, `page_faults_{major,minor}_median`
- `output_bytes_median`, `ok` (1 if all reps in cell exited cleanly)
- Per-stage medians (only those that have data appear):
  `total_run_time_s_median`, `sequence_1_io_s_median`,
  `seed_position_table_s_median`, `seed_hit_search_s_median`,
  `gapped_extension_s_median`, `chaining_s_median`, `output_s_median`,
  `ge_ydrop_one_sided_align_s_median` (the actual inner DP loop
  inside gapped extension — likely the GPU target),
  and a few other gapped-extension internal counters.

---

## 8. Stage-timing instrumentation

LASTZ's source already contains a `dbgTiming` framework (off by
default). Building with `-DdbgTiming -DdbgTimingGappedExtend` enables
accumulator-style timers around each algorithmic stage; on shutdown
LASTZ prints a per-stage report. Our `lastz_T` binary is exactly that.

**The three patches we made on top of upstream:**

1. The `dbg_timing_report` macro takes a `FILE*` so we can redirect the
   report.
2. The report is wrapped with `===STAGE_TIMING_BEGIN===` /
   `===STAGE_TIMING_END===` markers — deterministic for the parser to
   find inside stderr.
3. If `LASTZ_STAGE_REPORT=<path>` is set, the report is written to that
   file; otherwise it goes to stderr (upstream behavior).

`bench.py` sets `LASTZ_STAGE_REPORT=<run_dir>/runs/<cell>__rep<N>.stage.txt`
per rep, then parses it after the process exits. The parsed values
populate the `*_s_median` columns in `summary.csv`.

**Correctness:** instrumentation is pure bookkeeping (no control-flow
changes). `make test` produces byte-identical output between `lastz`
and `lastz_T` on the bundled smoke test.

**Overhead:** sub-millisecond on real workloads. To verify, run the
same suite with `--no-stage-timers` and compare wall medians.

---

## 9. Running on SLURM

**Don't run anything serious on the head node.** Even though lastz is
single-threaded, head-node CPU is shared with other users; numbers will
be noisy and you'll annoy people. Use `bench/sbatch/run_bench.sbatch`.

### Cluster config (current host)

```
sinfo:
  cpu  partition: 4 nodes (c0001..c0004), some drained
  gpu  partition: 4 nodes (g0001..g0004), 40 cores each, 4× Pascal GPU, 128 GB RAM
```

We submit to the `gpu` partition because that's where SegAlign will
need to run too — keeping baselines and GPU runs on the same hardware
removes a confounding variable.

### Submitting a suite

```bash
sbatch --export=ALL,SUITE=<suite>,REPS=<n>,RUN_NAME=<dir-name> \
       bench/sbatch/run_bench.sbatch
```

Examples:

```bash
# Bird-vs-bird chrZ, single rep, segalign_default variant.
sbatch --export=ALL,SUITE=bird_z_cross,REPS=1,RUN_NAME=bird_z_cross-2026-05-06 \
       bench/sbatch/run_bench.sbatch

# Same workload, isolate gapped-extension cost.
sbatch --export=ALL,SUITE=bird_z_cross_breakdown,REPS=1,RUN_NAME=bird_z_cross_breakdown-2026-05-06 \
       bench/sbatch/run_bench.sbatch

# Mammal pair (much heavier; bump time limit).
sbatch --time=06:00:00 --export=ALL,SUITE=,REPS=1,RUN_NAME=hg38_mm10-2026-05-06 \
       bench/sbatch/run_bench.sbatch   # add a suite for hg38_vs_mm10 first
```

The wrapper:
- Pins to the **first** CPU SLURM has actually given us
  (cgroup-aware — uses `taskset -pc $$` to discover allocated cpus,
  then takes the first one). You can override with `PIN_CPU=` in
  `--export`.
- Logs to `bench/sbatch/logs/lastz-bench.<jobid>.{out,err}`.
- Writes results to `bench/results/<RUN_NAME>/`.

### Monitoring

```bash
squeue -u $USER                                                       # what's running
tail -f bench/sbatch/logs/lastz-bench.<jobid>.err                     # bench.py progress
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State                  # post-mortem
scancel <jobid>                                                       # kill it
```

---

## 10. What we know so far

### Single-thread synthetic scaling (3-rep medians, default lastz)

From `bench/results/scaling-staged-v1/summary.csv`:

| Input size | Wall (s) | Peak RSS (MB) | seed_hit_search (s) | gapped_extension (s) | seed/total |
|---|---|---|---|---|---|
| 10 kbp × 10 kbp | 0.06 | 36 | 0.011 | 0.024 | 18% |
| 100 kbp × 100 kbp | 0.38 | 101 | 0.063 | 0.248 | 17% |
| 1 Mbp × 1 Mbp | 2.50 | 158 | 0.695 | 1.698 | 28% |
| 10 Mbp × 10 Mbp | 63.88 | 246 | **46.10** | 17.03 | **72%** |

**Two findings worth flagging:**

1. **Crossover.** At small inputs (≤1 Mbp) `gapped_extension` dominates;
   at 10 Mbp `seed_hit_search` overtakes and dominates (72%). This
   matches SegAlign's headline claim: at chromosome scale, the
   seed-and-filter stage is the bottleneck and the right thing to throw
   on a GPU.
2. **Super-linear scaling of seed_hit_search.** From 1 Mbp → 10 Mbp,
   input area grows 100×, total wall grows 25×, but
   `seed_hit_search` grows 66×. The likely cause is the diagonal-hash
   data structure (re-hashing + rehashed-bucket walks scale super-linearly
   with #seed-hits, which itself grows roughly with target_len ×
   query_len for random sequence).

### What we don't have yet

- **Anything at chromosome scale.** The bird-vs-bird run was just
  submitted (job 75052) at the time of writing — first chromosome-scale
  data point landing soon.
- **Any GPU comparison.** SegAlign is not yet built / run on this
  cluster; that's a follow-up.
- **Hardware counter data.** `perf` is unavailable on this kernel
  (paranoid=3, missing linux-tools for 5.4.0-216). Plan: source-level
  instrumentation with `__rdtsc`-style counters, or move to a host
  where we can run `perf stat`.

### Confirmed facts about LASTZ itself

- **Strictly single-threaded.** Zero pthread / OpenMP / MPI / fork
  code in `lastz/src/`. The way people scale lastz today is
  process-level sharding (split target by chromosome or by window, run
  N lastzes in parallel). SegAlign's wrapper does this with `eval &`.
  Implication: any GPU acceleration that wants to look credible should
  also report the N-core sharded baseline, not just 1-core.
- **Default scoring matrix is HOXD-derived** (4×4 nucleotide table
  hardcoded in `lastz.c`). `HOXD55` and `HOXD70` are *named alternate*
  matrices selectable via `--scores=`, but neither lastz nor SegAlign
  uses them by default.
- **SegAlign uses lastz defaults across the board.** Verified against
  their `src/main.cpp` lines 75–98:
  `--seed=12of19`, `--step=1`, `--xdrop=910`, `--hspthresh=3000`,
  `--ydrop=9430`, `--gappedthresh=hspthresh`. Their `--ydrop=9430` is
  literally `gap_open + 300 × gap_extend` evaluated against lastz's
  defaults. So `default` and `segalign_default` differ only in output
  format.

---

## 11. Open questions / next steps

In rough order:

1. **Land the chromosome-scale baseline.** Wait for job 75052
   (`bird_z_cross`) to finish; check `seed_hit_search` share. If >50%,
   it confirms the 10 Mbp synthetic finding scales up.
2. **Run `bird_z_cross_breakdown`** (default + nogapped on the same
   pair) to back out gapped-extension cost on real data, the same way
   we did for synthetic.
3. **Build + run SegAlign on the same chrZ pair** for a CPU-vs-GPU
   datapoint. Their docker image avoids the build complexity.
4. **Add an N-shard CPU baseline** to bench.py — launch N parallel
   `lastz` subprocesses on disjoint target slices (chromosome shards),
   produce a "CPU scaling curve" so the GPU comparison isn't 1-core vs
   1-GPU.
5. **Hardware counters.** Either fix `perf` (sysctl
   `kernel.perf_event_paranoid=2`) or instrument `seed_hit_search` and
   `ge_ydrop_one_sided_align` with `__rdtsc` counters for cycles + an
   estimate of memory-bandwidth utilization. We need to know whether
   the seed-hit hot loop is ALU- or memory-bound before designing a
   GPU kernel for it.
6. **Whole-genome runs.** Once chr-scale numbers look right, run
   hg38 vs mm10 whole-genome (overnight job). That's the SegAlign
   headline benchmark.

---

## 12. Gotchas, lessons learned, and "don't do that"

A list of things that wasted time during setup so they don't waste
time again:

- **Don't run benchmarks on the head node.** Use `sbatch`. Even
  single-thread numbers will be noisy under shared load.
- **`lastz` is single-threaded.** Allocating lots of CPUs for one
  benchmark process is wasteful; `--cpus-per-task=2` is enough (one for
  lastz, one for harness overhead). We use 4 to give the OS scheduler
  headroom.
- **`/usr/bin/time -v` keys contain colons.** "Elapsed (wall clock)
  time (h:mm:ss or m:ss): 0:01.23" — naive `split(":", 1)` mangles the
  key. `bench.py`'s `parse_time_v` uses `startswith` matching against a
  known key list. Don't refactor that without re-testing wall_s parsing.
- **Multi-sequence FASTAs need `[multiple]` action AND non-`lav`
  format.** lastz's default output format is `lav`, which can't
  represent multi-sequence inputs; `--format=maf` (or any non-lav) is
  required when `target_actions` includes `multiple`. Both
  `pseudopig_self` and `reads_vs_pseudopig` workloads hit this.
- **2bit subset syntax in lastz.** `lastz file.2bit[multiple,subset=chrZ]`
  does NOT work — lastz interprets `chrZ` as a literal filename. The
  correct way to pick one chromosome from a 2bit file is
  `lastz file.2bit/chrZ`. Documented in the lastz file-format reference
  but easy to miss.
- **UCSC `twoBitToFa` requires glibc 2.33+.** The host has glibc 2.27.
  Use `bench/twobit.py` instead (in-tree, pure Python, no deps).
- **NCBI accessions get bumped.** Initially picked an `XM_*` predicted
  variant for chimp ACTB and a stale fly Act5C accession; both
  returned HTTP 400 / wrong sequence. The fix is to use NCBI `esearch`
  + `esummary` to find a current curated `NM_*` accession, and have
  `fetch_real_data.py` carry a list `[primary, fallback, fallback2,...]`
  per target.
- **Header chars in FASTA length counts.** Original `seq_length_from_fasta`
  in `fetch_real_data.py` was including the `>NM_001101.5 ...` header
  characters in the length total, inflating reported lengths by ~80
  bytes. Fixed: skip lines starting with `>`.
- **`perf` is unusable here.** Kernel 5.4.0-216 with
  `kernel.perf_event_paranoid=3` and no `linux-tools` package. We rely
  on `/usr/bin/time -v` + in-binary stage timers. Don't waste time
  trying to build `perf` from source — use a different host or get
  sysadmin to drop the paranoid level.
- **Soft-masked vs unmasked input matters.** lastz's seeding behavior
  changes with `--unmask` / `--masking=N` flags. Our chrZ slices keep
  UCSC's soft-masking (lowercase = repeat); ~27% of galGal6.chrZ is
  lowercase. If you ever want hard-masked or unmasked variants, encode
  that as a workload flag, not a global edit.
- **`--ydrop=15000`, `HOXD55`, `--gappedthresh=3000` etc. are NOT
  SegAlign-specific.** SegAlign uses lastz defaults. If you see those
  values in old notes, ignore them — verified against the SegAlign
  source.

---

## Quick reference card

```bash
# Build
make -C lastz                       && make -C lastz build_lastz_timed

# Fetch data
python3 bench/fetch_real_data.py    # mRNAs (KB)
python3 bench/fetch_genomes.py      # genomes (GB, into scratch)

# Run smoke
python3 bench/bench.py --suite smoke --reps 3 --pin-cpu 0

# Run real (head-node OK only for tiny suites)
python3 bench/bench.py --suite actb_matrix --reps 3

# Run real (chromosome-scale → SLURM)
sbatch --export=ALL,SUITE=bird_z_cross,REPS=1,RUN_NAME=bird_z_cross-$(date +%F) \
       bench/sbatch/run_bench.sbatch

# Inspect
ls   bench/results/
cat  bench/results/<run>/summary.csv | column -ts,
cat  bench/results/<run>/runs/<wl>__<var>__rep0.stage.txt   # raw stage report
```

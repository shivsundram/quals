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
   - [5.4 Sub-chromosome slices (`slice_genome.py`)](#54-sub-chromosome-slices-slice_genomepy)
6. [Workload registry (`workloads.py`)](#6-workload-registry-workloadspy)
7. [Benchmark harness (`bench.py`)](#7-benchmark-harness-benchpy)
8. [Stage-timing instrumentation](#8-stage-timing-instrumentation)
9. [eu-stack PC sampling profiler (`profile_eustack.py`)](#9-eu-stack-pc-sampling-profiler-profile_eustackpy)
10. [Running on SLURM](#10-running-on-slurm)
11. [What we know so far](#11-what-we-know-so-far)
12. [Open questions / next steps](#12-open-questions--next-steps)
13. [Gotchas, lessons learned, and "don't do that"](#13-gotchas-lessons-learned-and-dont-do-that)

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
| `eu-stack` (from `elfutils`) | sampling profiler (§9) | `eu-stack --version` — `apt install elfutils` |
| `addr2line` (from binutils) | resolves PIE addrs to file:line for the profiler | `addr2line --version` |
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
    ├── fetch_genomes.py               # UCSC 2bit fetcher + per-chromosome FASTA extractor
    ├── twobit.py                      # in-tree .2bit reader (replaces UCSC twoBitToFa)
    ├── slice_genome.py                # sub-chromosome FASTA slicer (e.g. first 10 Mbp of chrZ)
    ├── profile_eustack.py             # eu-stack PC sampler — gives function- + file:line-level
    │                                  #   call-stack histograms for a single lastz run
    ├── data/                          # synthetic FASTAs (auto-created)
    ├── data_real/                     # NCBI fetches (auto-created)
    ├── data_genomes/                  # symlinks → /scratch2/.../2bit, slices, sub-slices
    ├── results/<run-name>/            # per-run output directory
    └── sbatch/
        ├── run_bench.sbatch           # SLURM wrapper — submits a suite to the gpu partition
        ├── profile_eustack.sbatch     # SLURM wrapper — submits one eu-stack profiling run
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
| `lastz/src/lastz_TS` | `make -C lastz build_lastz_substages` | **Substage-instrumented** (adds `-DdbgTimingSubstages` on top of `lastz_T`). Same wall-time output as `lastz_T` plus an rdtsc cycle breakdown of (i) the seed_hit_search hot loop (chain walk, processor callback, dedup-skip, x-drop, reporter) and (ii) the gapped-extension hot loop inside `ydrop_one_sided_align` (setup, row bounds, inner DP cells, row trailing, traceback, cleanup) — see §11 "Inside seed_hit_search at cycle granularity" and "Inside `ydrop_one_sided_align` at cycle granularity". Adds ~3-5 % overhead via per-iteration rdtsc on seed search and amortized <1 % overhead on gapped extension; use `lastz_T` for absolute wall-time, `lastz_TS` for substage breakdown. |

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
make -C lastz build_lastz_timed         # instrumented (lastz_T)
make -C lastz build_lastz_substages     # rdtsc substages (lastz_TS)
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

### 5.4 Sub-chromosome slices (`slice_genome.py`)

For profiler iteration we want runs that finish in tens of seconds, not
73 minutes. `slice_genome.py` cuts a contiguous range out of an existing
single-sequence FASTA (typically one of the per-chromosome files produced
by `fetch_genomes.py`), preserving N-blocks and soft-masking, and renames
the FASTA header so lastz sees a unique sequence id.

```bash
# First 10 Mbp of galGal6.chrZ
python3 bench/slice_genome.py \
    bench/data_genomes/galGal6.chrZ.fa \
    bench/data_genomes/galGal6.chrZ_0_10mb.fa \
    --start 0 --length 10_000_000 \
    --rename "galGal6.chrZ:0-10mb"

# First 10 Mbp of taeGut2.chrZ
python3 bench/slice_genome.py \
    bench/data_genomes/taeGut2.chrZ.fa \
    bench/data_genomes/taeGut2.chrZ_0_10mb.fa \
    --start 0 --length 10_000_000 \
    --rename "taeGut2.chrZ:0-10mb"
```

The script prints a small stats block (ACGT / acgt / N counts) so you can
sanity-check the masking density of your slice — e.g. `galGal6.chrZ:0-10mb`
is 21.0% lowercase (vs 27% for the whole chrZ — the 5' tip is unusually
clean) and `taeGut2.chrZ:0-10mb` is only 11.7% lowercase. Worth knowing
when comparing slices to whole-chr numbers.

These slices live under `/scratch2/shiv1/lastz-bench-data/slices/` and are
symlinked into `bench/data_genomes/` exactly like the full-chr files.

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
real.galGal6_vs_taeGut2_chrZ_10mb      first 10 Mbp of each chrZ    ← profiler iteration target
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
bird_z_cross_10mb           first 10 Mbp of each chrZ, segalign_default     ← profiler iteration
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

**Overhead:** sub-millisecond on real workloads in aggregate. The eu-stack
profile (§9) shows ~3-4% of wall time inside `__vdso_gettimeofday`, so
the timers themselves are visible but not dominant. For absolute
wall-time numbers, sanity-check with `--no-stage-timers` (uses the
non-instrumented `lastz` binary).

`lastz_T` is **also built with `-g`** (DWARF debug info) so eu-stack and
addr2line can resolve hot addresses to file:line. `-g` is free at `-O3`
— it only adds `.debug_*` ELF sections, which the loader doesn't touch
and which don't affect `.text` size or branch layout.

---

## 9. eu-stack PC sampling profiler (`profile_eustack.py`)

We can't use `perf` on this cluster (kernel 5.4 +
`kernel.perf_event_paranoid=3`, no `linux-tools` package). The next
best thing is **`eu-stack` from elfutils**: a stack walker that, when
called repeatedly, becomes a poor-man's sampling profiler. Combined
with `addr2line` against `lastz_T` (built with `-g`), we get function-
**and** file:line-level breakdown of where the program counter is
spending its time.

`bench/profile_eustack.py` wraps the whole thing:

1. Forks `lastz_T` as a child (so we own its ptrace permission under
   Ubuntu's default `ptrace_scope=1`) and has the child call
   `prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY)` via `preexec_fn`, so
   sibling processes (the `eu-stack` invocations) can attach.
2. Captures the PIE load base from `/proc/<pid>/maps` while the target
   is alive.
3. Loops `eu-stack -p <pid>` at the requested rate, saving each raw
   sample.
4. After the target exits, runs `addr2line` on every unique leaf PC
   (with the load-base subtracted, since the binary is PIE), and
   aggregates into:
   - `folded.txt` — Brendan-Gregg-style stack folding
     (`outer;...;inner <count>`), ready to feed to `flamegraph.pl`.
   - `top_funcs.txt` — flat function-level histogram.
   - `top_lines.txt` — leaf function + file:line + count, line-resolved
     where DWARF info exists, falling back to hex address otherwise.
   - `manifest.json` — what we ran, when, sample count, load base.

### CLI

```bash
python3 bench/profile_eustack.py \
    --lastz   lastz/src/lastz_T \
    --target  bench/data_genomes/galGal6.chrZ_0_10mb.fa \
    --query   bench/data_genomes/taeGut2.chrZ_0_10mb.fa \
    --out-dir bench/results/eustack-bird10mb-001 \
    --rate-hz 100 \
    --extra-arg=--format=maf-
```

Notes:

- Use `--extra-arg=<value>` (with `=`) for any lastz arg that starts
  with `--`. argparse won't accept it as a separate token.
- `--rate-hz 100` (default) gives ~5 ms between samples, which is well
  below the per-sample eu-stack cost (~1-3 ms on this host). At higher
  rates the sampler can't keep up and effective rate drops.
- `--max-samples N` caps the number of samples (useful for very long
  runs where you don't want hundreds of MB of raw stacks).
- Default keeps `raw/sample-NNNNNN.txt` for re-aggregation; use
  `--no-keep-raw` to discard them after folding.

### Sizing

| Workload | Wall (s) | Samples @ 100 Hz | Stat resolution |
|---|---|---|---|
| synthetic 1 Mbp × 1 Mbp | ~3 | ~300 | OK for top-3 functions |
| 10 Mbp × 10 Mbp slice | ~60-90 | ~6000-9000 | Good (1% resolution) |
| Full bird-Z chr-cross | ~73 min | ~440k @ 100 Hz | Overkill — use lower rate or `--max-samples` |

### Smoke result (synthetic 1 Mbp × 1 Mbp, default lastz, 312 samples)

```
ydrop_one_sided_align.part.0          57.7%   (gapped extension inner DP loop)
find_table_matches                    18.3%   (seed lookup; hot at seeds.c:1306)
xdrop_extend_seed_hit                  9.3%   (ungapped extension after a seed hit)
__vdso_gettimeofday                    3.8%   (stage-timer overhead — see §8)
gapped_extend                          3.5%
process_for_simple_hit                 2.2%
add_word                               1.9%   (seed table build)
```

Use this as the synthetic-DNA baseline; the bird-Z 10 Mbp profile is
the real-DNA equivalent and is the next thing to capture (see §12).

### Submitting a profile run via SLURM

For anything beyond a few seconds of wall time, run on a compute node:

```bash
sbatch --export=ALL,\
    TARGET=bench/data_genomes/galGal6.chrZ_0_10mb.fa,\
    QUERY=bench/data_genomes/taeGut2.chrZ_0_10mb.fa,\
    OUT_DIR=bench/results/eustack-bird10mb-001 \
    bench/sbatch/profile_eustack.sbatch
```

The wrapper accepts `RATE_HZ`, `EXTRA_ARG`, `LASTZ`, `PIN_CPU` overrides
(see comments in `bench/sbatch/profile_eustack.sbatch`). It pins both
the python sampler and the lastz child to a single CPU via `taskset`,
just like the bench harness wrapper.

### Visualizing the results

```bash
# Quick text summary
cat bench/results/eustack-bird10mb-001/top_funcs.txt | head -15
cat bench/results/eustack-bird10mb-001/top_lines.txt | head -25
head -10 bench/results/eustack-bird10mb-001/folded.txt

# Flamegraph (if you have Brendan Gregg's flamegraph.pl in PATH)
flamegraph.pl < bench/results/eustack-bird10mb-001/folded.txt > flame.svg

# Re-aggregate from raw stacks (e.g. to filter out specific frames):
python3 -c '
import sys, collections
counts = collections.Counter()
for path in sys.argv[1:]:
    with open(path) as f:
        frames = [line.split()[2] for line in f if line.lstrip().startswith("#")]
    if frames:
        counts[";".join(reversed(frames))] += 1
for k, v in counts.most_common():
    print(v, k)
' bench/results/eustack-bird10mb-001/raw/sample-*.txt
```

### Known limitations

- **No GPU profiling** — eu-stack is CPU-only.
- **PIE addresses are randomized per run.** We capture the load base
  from `/proc/<pid>/maps` while the target is alive and use it for
  addr2line; if `/proc/<pid>/maps` can't be read in time (very short
  runs, <50 ms), line resolution silently falls back to hex addresses.
- **Inlined functions get attributed to their callsite, not their
  definition.** GCC at `-O3` partial-inlines `ydrop_one_sided_align`
  into multiple seed_search.c callsites, so the "file:line" for
  `ydrop_one_sided_align.part.0` shows up as `seed_search.c:399`
  (the call site) rather than `gapped_extend.c:NNN`. eu-stack `-i`
  (show inlined frames) helps but isn't perfect; we don't enable it by
  default to keep the folded output simpler.
- **Effective sample rate ≤ requested rate.** eu-stack costs 1-3 ms
  per call (process spawn + ptrace attach + walk + detach). At 200 Hz
  requested we typically see 75-100 Hz effective. Bump SLURM cpus if
  you want sustained 200+ Hz.
- **ptrace_scope=1 is fine; ptrace_scope=2 or 3 is not.** If you ever
  see `dwfl_thread_getframes: Operation not permitted` after our
  PR_SET_PTRACER opt-in, the host has been locked down further; ask
  the sysadmin to lower `kernel.yama.ptrace_scope`.

---

## 10. Running on SLURM

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

# eu-stack profile of the bird-Z 10 Mbp slice.
sbatch --export=ALL,\
    TARGET=bench/data_genomes/galGal6.chrZ_0_10mb.fa,\
    QUERY=bench/data_genomes/taeGut2.chrZ_0_10mb.fa,\
    OUT_DIR=bench/results/eustack-bird10mb-001 \
    bench/sbatch/profile_eustack.sbatch
```

The bench wrapper (`run_bench.sbatch`):
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

## 11. What we know so far

### Single-thread synthetic scaling (3-rep medians, default lastz)

From `bench/results/scaling-staged-v1/summary.csv`:

| Input size | Wall (s) | Peak RSS (MB) | seed_hit_search (s) | gapped_extension (s) | seed/total |
|---|---|---|---|---|---|
| 10 kbp × 10 kbp | 0.06 | 36 | 0.011 | 0.024 | 18% |
| 100 kbp × 100 kbp | 0.38 | 101 | 0.063 | 0.248 | 17% |
| 1 Mbp × 1 Mbp | 2.50 | 158 | 0.695 | 1.698 | 28% |
| 10 Mbp × 10 Mbp | 63.88 | 246 | **46.10** | 17.03 | **72%** |

### Real chromosome-scale: two contrasting profiles

Both runs use the `segalign_default` variant (lastz defaults +
`--format=maf-`), single thread.

**galGal6.chrZ × taeGut2.chrZ** — bird-vs-bird cross, 82 Mbp × 73 Mbp,
real evolutionary divergence (~100 Mya). On g0002 (GPU node, pinned).
From `bench/results/bird_z_cross-2026-05-06/.../__-1.stage.txt`:

| Stage | Wall (s) | % of total |
|---|---|---|
| **Total** | **4366** (73 min) | 100% |
| seed_hit_search | **3995** | **91.5%** |
| gapped_extension | 363 | 8.3% |
| I/O + index + output | 7 | 0.16% |

**galGal6.chrZ × galGal6.chrZ** — same chromosome aligned to itself,
82 Mbp × 82 Mbp. On sapling2 (head node, pinned). From
`bench/results/bird_z_self-2026-05-06/.../__0.stage.txt`:

| Stage | Wall (s) | % of total |
|---|---|---|
| **Total** | **14346** (239 min, ~4 h) | 100% |
| seed_hit_search | 2532 | 17.6% |
| **gapped_extension** | **11800** | **82.2%** |
| I/O + index + output | 14 | 0.10% |
| Output produced | **1.82 GB** of MAF | (vs 61 MB for cross) |

In *both* runs, inside gapped_extension, `ydrop_one_sided_align` (the
inner banded-SW DP loop) accounts for ~96% of the gapped time:
349 s of 363 s (cross) and 5336 s of 11800 s (self).

**The seed/gapped ratio inverts based on alignment shape.** Self-vs-self
produces a giant trivial diagonal (the same chromosome aligning to
itself perfectly along its full length) — that single alignment alone
fills most of the gapped-extension budget and produces 30× more output
bytes than the cross run. Cross-species at moderate divergence produces
many more spurious seed hits that get ungapped-extended and rejected,
filling the seed_hit_search budget.

For the SegAlign comparison the **cross profile is the relevant one** —
nobody benchmarks whole-genome self-self alignment. But the self
profile is a useful reminder that "what's the bottleneck" depends
strongly on alignment density, not just input size.

### Inside seed_hit_search: the bird-Z 10 Mbp eu-stack profile

The chr-scale stage report tells us `seed_hit_search` is 91.5% of wall
time on bird-Z cross — but not which **part** of `seed_hit_search`. To
break that out, we ran `bench/profile_eustack.py` on the 10 Mbp slice
(`real.galGal6_vs_taeGut2_chrZ_10mb`, ~56 s wall, 3139 samples at
100 Hz). Output in `bench/results/eustack-bird10mb-001/`.

**Top exclusive frames (where the IP literally was at sample time):**

| Function | Exclusive % | What it does |
|---|---:|---|
| `xdrop_extend_seed_hit` | **48.6 %** | ungapped x-drop extension from a seed hit |
| `find_table_matches` | **34.7 %** | chain walk through `pt->prev[]` for one query word |
| `ydrop_one_sided_align` | 11.0 % | gapped-extension inner DP loop |
| `process_for_simple_hit` | 3.9 % | diagonal-dedup gatekeeper that calls x-drop |
| everything else | 1.8 % | seeding, I/O, parsing, stage-timer overhead |

**Inclusive call tree** (the same data viewed as cumulative time per
call subtree — what the flamegraph shows). Inclusive % includes the
function itself plus everything it calls; exclusive % is just the body:

```
seed_hit_search loop body                    (~0 % exclusive — all in callees)
└─ for each surviving seed position in target chain:                 ←─ 87 % inclusive
     │                                                                   here
     ├─ chain walk: load pt->prev[s-1]      ────┐
     │   advance s, load DNA, check word      ├── 35 % EXCL
     │   = find_table_matches body            └─  in find_table_matches
     │
     └─ process_for_simple_hit callback        ←── 52 % inclusive
          ├─ diagonal dedup check                  4 % EXCL
          │   diagEnd[h] = compare endpoints       in process_for_simple_hit
          │
          └─ xdrop_extend_seed_hit              ←── 49 % EXCL = leaf
              ├─ left  extension scoring loop
              ├─ right extension scoring loop
              └─ update diagEnd[h] for next time
```

Read the flamegraph as: every box's width = **inclusive** time. The
strip of a parent box that has nothing stacked on top of it = its
**exclusive** time. So the bottom red bar (`find_table_matches`) at
~87 % width has a ~52 %-wide `process_for_simple_hit` stacked on it,
leaving a ~35 % strip with no children — that strip is the chain walk
loop body itself.

**What's in the 35 % "exclusive find_table_matches"?** That's the
inner loop body at `lastz/src/seed_search.c:832-872`:

```
for (pos = pt->last[packed2]; pos != noPreviousPos; pos = pt->prev[pos]) {
    pos1 = adjStart + step*pos;          // arithmetic — cheap
    if ((selfCompare) && ...) continue;  // false on cross-species, ~free
    if ((sameStrand)  && ...) continue;  // bandWidth==0 default, short-circuits
    seed_search_count_stat(rawSeedHits); // no-op when stats off
    basesHit += (*processor)(...);       // → process_for_simple_hit
}
```

The 35 % is dominated by **one thing**: the pointer-chase load
`pos = pt->prev[pos]`. Each iteration's address depends on the
previous iteration's load result, so the CPU can't speculate ahead and
can't usefully prefetch — `pt->prev[]` is a flat `u32` array indexed
by an essentially-random target position (whatever happened to hash
to the same packed seed word as our query word). On real DNA with
soft-masked repeats, those chains can be long (a popular 12-mer can
appear thousands of times across a chromosome), and walking each one
is a serial L2/L3 latency ladder.

Secondary contributors inside the 35 %:
- **Indirect call dispatch** through `*processor` — GCC can't inline
  `process_for_simple_hit` because the call goes through a function
  pointer (`hp->hitProc`) set up at workload type dispatch time. So we
  pay a `call/ret` per seed hit even when dedup short-circuits.
- **A multiply + add** for the position translation
  (`pos1 = adjStart + step*pos`).
- **Loop overhead** (compare against `noPreviousPos`, branch).

We don't have a finer cycle breakdown yet — that's the rdtsc-substages
patch's job (next pending todo). The expected outcome is that
~70-90 % of the 35 % is the `pt->prev[]` load latency, with the rest
in indirect-call dispatch and arithmetic. If true, **the chain walk
is memory-latency-bound** and the right GPU port pattern is what
SegAlign does: flatten the `pt->last`/`pt->prev` linked list into a
dense per-bucket array so the chain walk becomes a contiguous stream
load, friendly to wide vector loads + bandwidth.

**Findings from the eu-stack profile that the stage report missed:**

1. The dominant single function is `xdrop_extend_seed_hit` at 48.6 %,
   not `find_table_matches`. Roughly half the entire program at chr
   scale is the ungapped x-drop scoring loop — exactly the kernel
   SegAlign accelerates first on the GPU.
2. The seed-lookup-vs-ungapped-extend ratio inverts between scales:
   1 Mbp synth = 1.97 : 1 (lookup-heavy), 10 Mbp real = 0.71 : 1
   (extend-heavy). At chr scale, more seed hits **survive** the dedup
   filter and reach the x-drop call.
3. Diagonal dedup is **invisible** in the profile (the dedup check
   itself is hidden inside the 4 % exclusive `process_for_simple_hit`
   strip). That's a feature, not a bug — it means the cache-line
   `diagEnd[h]` lookup is essentially free, and it's filtering enough
   redundant extensions that the savings dwarf its own cost.
4. Stage-timer overhead is amortized away at chr scale: 3.8 % at 1 Mbp
   synth → 0.1 % at 10 Mbp real. The instrumented `lastz_T` binary is
   safe to use as the production profiling target — no need to also
   re-run with `--no-stage-timers`.
5. Only **14 distinct call stacks** across 3139 samples. lastz spends
   essentially all of its time on a handful of code paths, which is
   excellent news for porting: cover 2-3 functions and you've covered
   ~90 % of the runtime.

### Inside seed_hit_search at cycle granularity (rdtsc substages)

The eu-stack profile gives function-level *time* shares; it cannot tell us
whether the chain walk is memory-latency-bound or compute-bound, or how
much of the x-drop bar is failed extension attempts vs. successful HSPs.
For that we built `lastz_TS` (`make -C lastz build_lastz_substages`),
which adds `__rdtsc()` cycle accumulators around four hot blocks in
`lastz/src/seed_search.c`:

- `ftm.chain_walk` — chain-walk loop body in `find_table_matches` minus
  the time inside the processor callback. This is the pure pointer
  chase through `pt->prev[pos]`.
- `ftm.processor_callback` — total cycles in `(*processor)(...)`, which
  for the `segalign_default` workload always dispatches to
  `process_for_simple_hit`.
- `pfs.dedup_short_circuit` — calls that returned via the `diagEnd[hDiag]`
  dedup check, plus the cycles spent in that path.
- `xdrop_extend_seed_hit` — calls and cycles inside the ungapped
  extension, plus a count of how many returned `noScore` (i.e. failed
  to find an HSP).
- `reporter_callback` — calls and cycles in `(*info->hp.reporter)(...)`,
  the HSP-emit dispatch.

Run the same workload through `lastz_TS`; the report appears inside the
existing `===STAGE_TIMING_BEGIN===` / `===STAGE_TIMING_END===` block, so
no harness changes are needed. Raw output of both runs in
`bench/results/rdtsc-substages-bird10mb-001/`.

**Synthetic 1 Mbp × 1 Mbp** (3.0 s wall, 26 M find_table_matches calls,
2.1 M chain iters):

| Substage | Total cyc | Calls/iters | Cyc/call | Notes |
|---|---:|---:|---:|---|
| `ftm.chain_walk` | 186 M | 2.10 M | **88.9** | ~30 ns @ 3 GHz, mostly L2/L3 (table fits) |
| `ftm.processor_callback` | 993 M | 2.10 M | 473.6 | dominated by x-drop subtree |
| `pfs.dedup_short_circuit` | 16 M | 542 k (26 % of pfs) | 30.3 | one cache-line lookup + branch |
| `xdrop_extend_seed_hit` | 820 M | 1.55 M | 527.7 | of which 99.5 % return noScore |
| `reporter_callback` | 1.4 M | 6,744 | 210.8 | rare — only successful HSPs |

**Real bird-Z 10 Mbp × 10 Mbp** (62.3 s wall, 226 M find_table_matches
calls, 165 M chain iters):

| Substage | Total cyc | Calls/iters | Cyc/call | Notes |
|---|---:|---:|---:|---|
| `ftm.chain_walk` | 20.6 G | 165 M | **125.3** | **+40 % vs 1 Mbp synth** — DRAM territory |
| `ftm.processor_callback` | 81.8 G | 165 M | 496.1 | per-call cost flat across scales |
| `pfs.dedup_short_circuit` | 67 M | 1.85 M (**1.1 %** of pfs) | 36.0 | barely fires on real DNA! |
| `xdrop_extend_seed_hit` | 69.2 G | 163 M | 424.8 | **99.999 %** return noScore (1665 hits / 163 M attempts) |
| `reporter_callback` | 0.28 M | 1665 | 170.3 | the HSPs the eu-stack flame ends in |

**Three findings the rdtsc data reveals that eu-stack alone could not:**

1. **Chain walk is memory-latency-bound.** Cyc-per-iter grew 89 → 125
   (+40 %) going from 1 Mbp synth to 10 Mbp real. At ~3 GHz that's
   30 ns → 42 ns per pointer chase, sitting between L3 hit latency
   (~12 ns) and DRAM (~80 ns). At 10 Mbp the `prev[]` array is ~60 MB
   (well past L3) and access patterns are random, so most loads miss
   to DRAM. **This is the quantitative confirmation that switching to
   the dense-array layout (counting-sort buckets) is worth doing on
   the CPU before going GPU**: same algorithm, but reads become
   contiguous and the hardware prefetcher can hide most of the
   latency. Expected savings: most of the 7 s currently in chain walk.
2. **Diagonal dedup barely fires on real cross-species DNA.** The
   `diagEnd[hDiag]` short-circuit catches **26 %** of hits on
   synthetic 1 Mbp but only **1.1 %** on real bird-Z 10 Mbp. Why:
   soft-masking already suppresses the repetitive seed hits that
   would create redundant diagonals. Implication: **SegAlign throwing
   away dedup on the GPU costs them ~1 % of x-drop work on real
   cross-species, not the 26 % the synthetic profile suggested.** A
   GPU port doesn't need to faithfully reproduce dedup; it can rely on
   masking + sheer parallelism instead.
3. **99.999 % of x-drop calls return noScore — 162.97 M failures vs.
   1665 successes.** Of the 23 s spent in `xdrop_extend_seed_hit`,
   almost all of it is paying 425 cycles per call to discover "no,
   this seed hit does not extend to a real HSP". The CPU is doing
   163 M serially-dependent failed-extension scoring attempts. **This
   is the strongest argument for GPU x-drop acceleration**: 163 M
   nearly-independent, identical-shape compute units at 0.001 %
   success rate. Embarrassingly parallel; the CPU pays 425 cycles
   serially for every "no" answer where a GPU pays the same on every
   warp lane in parallel.

**Reconciliation check.** ftm.processor_cyc (81.8 G) − pfs subtree sum
(69.3 G = dedup + xdrop + reporter) = 12.5 G cyc gap = 15 % of the
processor budget. That's the body of `process_for_simple_hit` itself —
the `hashedDiag` math, the `gfExtend` dispatch, and call/ret overhead.
Same ~15 % gap on synthetic; the bookkeeping is consistent across
workloads.

#### Seed-weight sweep — sensitivity vs. speed on the same workload

Quick A/B/C with three seed shapes on bird-Z 10 Mbp (full numbers and
analysis in `bench/results/seed-sweep-bird10mb/SUMMARY.md`):

| Variant | Wall | Chain iters | HSPs | Sens. vs default |
|---|---:|---:|---:|---:|
| `--seed=12of19` (default) | 62.3 s | 164.8 M | 510 | 100 % |
| `--seed=14of22` (heavy spaced) | **17.1 s** (3.6×) | 13.2 M | 490 | **96.1 %** |
| `--seed=match12` (contiguous) | **10.4 s** (6.0×) | 18.7 M | 436 | **85.5 %** |

What the cycle-level data adds:

- The **chain-walk cyc/iter is identical** for `12of19` and `match12`
  (125 cyc/iter both) — the DRAM-latency cost is set by the table
  layout, not the seed shape. Heavier seeds shrink the absolute hit
  count but don't fix the memory pattern; the dense-array layout fix
  is orthogonal to seed weight.
- **Dedup short-circuits 27 % of `match12` calls but only 1 % of
  `12of19` calls.** Diagonal dedup was designed for contiguous seeds;
  spaced seeds' don't-care positions scatter repeats across many
  buckets and do dedup's job before lookup. SegAlign throwing dedup
  away costs even less than the 1 % we measured for the default seed.
- The **99.99 % x-drop failure rate is structural**, not seed-dependent.
  All three variants land at >99.98 % `noScore`. Heavier seeds reduce
  the absolute count of failed extensions but don't change the success
  rate. The GPU's "embarrassingly parallel failure" win applies at any
  seed weight.

Practical takeaway: `14of22` is a near-Pareto win for cross-species
alignment at vertebrate-scale divergence (~3.6× speedup, ~4 %
sensitivity loss). SegAlign keeping `12of19` as the safe default makes
sense for "we want to handle anything from primates to fish" but
heavier seeds are the right pick for closer species. **Seed weight is
a sensitivity dial, not an algorithmic fix** — the dense-array layout
+ GPU x-drop wins remain the structural ones.

**Reproduce.** Prereqs: `lastz_TS` built (`make -C lastz
build_lastz_substages`), and the bird-Z 10 Mbp slices already on disk
(`bench/data_genomes/galGal6.chrZ_0_10mb.fa` and
`taeGut2.chrZ_0_10mb.fa` — produced earlier; see §5.4).

```bash
mkdir -p bench/results/seed-sweep-bird10mb
cd bench/results/seed-sweep-bird10mb

TARGET=$(realpath ../../../bench/data_genomes/galGal6.chrZ_0_10mb.fa)
QUERY=$(realpath ../../../bench/data_genomes/taeGut2.chrZ_0_10mb.fa)
LASTZ=$(realpath ../../../lastz/src/lastz_TS)

run_one() {
    local tag="$1"; shift
    /usr/bin/time -v -o ${tag}.time.txt \
        env LASTZ_STAGE_REPORT=${tag}.stage.txt \
            taskset -c 0 \
            "$LASTZ" "$TARGET" "$QUERY" "$@" --format=maf- \
            > ${tag}.maf 2>${tag}.err
    local hsp=$(grep -c '^a score=' ${tag}.maf)
    local wall=$(grep 'Elapsed (wall' ${tag}.time.txt | awk '{print $NF}')
    local seedhit=$(grep '^seed hit search:' ${tag}.stage.txt | awk '{print $NF}')
    echo "$tag: wall=$wall  seed_hit=${seedhit}s  HSPs=$hsp"
}

run_one default_12of19                       # bare command, default --seed=12of19
run_one heavy_14of22       --seed=14of22
run_one contig_match12     --seed=match12
```

Each variant takes ~1 minute (default), ~17 s (heavy), or ~10 s
(contig). Total ~90 s on a head-node CPU 0. The rdtsc substage block
appears between `===STAGE_TIMING_BEGIN===` / `===STAGE_TIMING_END===`
markers in `<tag>.stage.txt`; pull it out with:

```bash
sed -n '/--- rdtsc substage breakdown/,/===STAGE_TIMING_END/p' \
    default_12of19.stage.txt
```

#### Within-species regime change: hg38 vs T2T-CHM13 chr1:50M-60M

To pressure-test the seed-tuning intuition at the *opposite* end of
the divergence spectrum, we ran the same 10 Mbp × 10 Mbp matrix using
two human assemblies (~99.5–99.9 % identity in single-copy regions) at
`chr1:50,000,000-60,000,000`. Three variants, same `lastz_TS` binary,
single-thread CPU 0. Full numbers in
`bench/results/seed-sweep-hg-vs-chm13/SUMMARY.md`.

| Variant | Wall | seed_hit | gapped_ext | HSPs | bp aligned |
|---|---:|---:|---:|---:|---:|
| `--seed=12of19` (cross-species default) | 47.6 s | 18.2 s | 28.5 s | 2474 | 12.44 Mbp |
| `--seed=match15` (bigger contig seed) | 26.8 s | **0.83 s** | 24.9 s | 1798 | 12.12 Mbp |
| `--seed=match15 --step=20 --notransition` | 21.5 s | **0.37 s** | 20.5 s | 943 | 11.33 Mbp |

Three findings worth pulling out:

**1. The bottleneck *flips* at high identity.** On bird-Z 10 Mbp cross-
species (~75 % identity), `seed_hit_search` was 88 % of total wall
time. On hg38-vs-CHM13 with cross-species defaults it's already only
38 %; with `--seed=match15 --step=20` it drops to 2 %. **Gapped
extension goes from 10 % to 95 % of runtime.** Within-species and
cross-species are not the same problem. SegAlign's GPU port of *both*
seed_hit_search and gapped extension is essential — accelerating only
the seed lookup would leave ~95 % of the runtime on the CPU for
within-species workflows.

**2. HSP count is a misleading sensitivity metric.** Going from
`12of19` to `match15+step=20+notransition` drops HSP count by 62 %
(2474 → 943) but **total bp aligned drops only 9 %** (12.44 → 11.33
Mbp), and **mean HSP length more than doubles** (5,029 → 12,014 bp).
Maximum HSP length is identical (~448 kbp) across all three. The
"missing" HSPs aren't missing homology — they're short fragments that
get reported as a single longer HSP when the seed grid is sparser.
**Use total-bp-aligned, not HSP count, to compare sensitivity across
seed configurations.**

**3. seed_hit_search wall savings don't translate 1:1 to total wall
speedup.** seed_hit_search drops 49× (18.2 s → 0.37 s) but total wall
only drops 2.2× (47.6 s → 21.5 s). Why: gapped extension is the new
bottleneck and it shrinks only ~30 % (28.5 s → 20.5 s) — set by the
amount of homology found, not by the seed machinery. **To get bigger
wall wins on within-species data, the GPU must also accelerate gapped
extension.** This is exactly what SegAlign does (its `processSegment`
kernel does banded x-drop on GPU); the seed-only fix is necessary but
not sufficient at this regime.

Combined cross-species + within-species picture for the GPU port:

| Regime | Bird-Z cross-species | hg-vs-CHM13 within-species |
|---|---:|---:|
| Identity | ~75 % | ~99.5 % |
| seed_hit_search share | 88 % | 38 % → 2 % with tuning |
| gapped_extension share | 10 % | 60 % → 95 % with tuning |
| Where does GPU help most? | Chain walk + x-drop kernels | Gapped extension kernel |

**The takeaway is structural**: depending on the use case, "where does
the time go" shifts dramatically. A serious GPU-accelerated lastz needs
to attack both bottlenecks because the same tool serves both regimes.

Caveat on rdtsc coverage: `match15` has bit-weight 30, which exceeds
the 24-bit max packed-index size, so lastz takes the "overweight seed"
path through `find_table_matches_resolve()` — a different function
than the `find_table_matches()` we instrumented for the bird-Z
sweep. The chain-walk substage timer reads 0 on these runs as a
result (the work happens in the un-instrumented function); xdrop and
reporter timers are unaffected. Adding the same `__rdtsc()` block to
`find_table_matches_resolve` is a small follow-up if we want full
substage coverage on overweight seeds.

**Reproduce.** Prereqs: `lastz_TS` built (`make -C lastz
build_lastz_substages`); hs1 (T2T-CHM13) fetched and chr1 extracted
from both human assemblies. The data fetch is one-time and ~30 s for
the download itself plus ~20 s for the chr1 extraction:

```bash
# 1) Fetch T2T-CHM13 (UCSC code "hs1") — ~775 MB 2bit, ~30 s download.
#    Skips if already cached at /scratch2/shiv1/lastz-bench-data/2bit/.
python3 bench/fetch_genomes.py --only hs1

# 2) Extract chr1 from hg38 (only chr19 was extracted by default) and
#    symlink alongside hs1.chr1.fa which fetch_genomes.py made above.
python3 -m bench.twobit extract \
    /scratch2/shiv1/lastz-bench-data/2bit/hg38.2bit chr1 \
    /scratch2/shiv1/lastz-bench-data/slices/hg38.chr1.fa
ln -sf /scratch2/shiv1/lastz-bench-data/slices/hg38.chr1.fa \
    bench/data_genomes/hg38.chr1.fa

# 3) Slice 10 Mbp at chr1:50M-60M (gene-rich p-arm, no N-blocks here).
python3 bench/slice_genome.py \
    bench/data_genomes/hg38.chr1.fa \
    bench/data_genomes/hg38.chr1_50_60mb.fa \
    --start 50000000 --length 10000000 --rename hg38_chr1_50M_60M
python3 bench/slice_genome.py \
    bench/data_genomes/hs1.chr1.fa \
    bench/data_genomes/hs1.chr1_50_60mb.fa \
    --start 50000000 --length 10000000 --rename hs1_chr1_50M_60M
```

Then the sweep itself (~95 s wall, single-thread CPU 0):

```bash
mkdir -p bench/results/seed-sweep-hg-vs-chm13
cd bench/results/seed-sweep-hg-vs-chm13

TARGET=$(realpath ../../../bench/data_genomes/hg38.chr1_50_60mb.fa)
QUERY=$(realpath  ../../../bench/data_genomes/hs1.chr1_50_60mb.fa)
LASTZ=$(realpath  ../../../lastz/src/lastz_TS)

run_one() {
    local tag="$1"; shift
    /usr/bin/time -v -o ${tag}.time.txt \
        env LASTZ_STAGE_REPORT=${tag}.stage.txt \
            taskset -c 0 \
            "$LASTZ" "$TARGET" "$QUERY" "$@" --format=maf- \
            > ${tag}.maf 2>${tag}.err
}

run_one cross_species_default                                              # default --seed=12of19
run_one bigger_seed_match15        --seed=match15
run_one hg_recommended             --seed=match15 --step=20 --notransition
```

The sensitivity check (HSPs vs total bp aligned) is a one-liner over
the resulting MAF blocks:

```bash. 
for f in *.maf; do
    awk -v tag=${f%.maf} '
        /^a score=/ { getline; split($0, a, " "); sum += a[4]; n++ }
        END { printf "%-25s HSPs=%d total_bp=%d mean=%d\n",
                     tag, n, sum, n ? sum/n : 0 }' "$f"
done
```

Use `total_bp` as the sensitivity metric, not HSP count — see finding
#2 above for why.

**Cycle accounting against wall time.** Bird-Z 10 Mbp at 62.3 s:
- Chain walk: 20.6 G cyc → ~7 s (assuming 3 GHz; could be less at
  turbo throttling)
- Processor callback total: 81.8 G cyc → ~27 s
- Outside seed_hit_search (gapped extension etc.): 6.2 s by stage
  timer
- Total ≈ 40-50 s of accounted CPU time vs. 62 s wall. The remainder
  is dead reckoning between rdtsc (true cycles) and wall time (turbo
  + memory stalls counted as wall but not as retired cycles), and is
  consistent with what we'd expect from a CPU running closer to its
  base clock under heavy DRAM traffic.

#### Inside `ydrop_one_sided_align` at cycle granularity (rdtsc substages)

In four of five regimes, gapped extension is 83-99 % of wall (sparse
mouse-human, dense bird-Z, dense mouse-rat, dense human-rhesus, dense
within-species — see "Step-axis Pareto sweeps" below). The existing
`-DdbgTimingGappedExtend` instrumentation tells us that virtually all
of that time is in `ydrop_one_sided_align` — but not where inside it.
We extended the same `lastz_TS` build with rdtsc accumulators around
six substages of `ydrop_one_sided_align`:

| substage | what it covers | granularity |
|---|---|---|
| `yda.setup`         | scoring extract, L/R bound init, first DP row     | per call |
| `yda.row_bounds`    | `update_LR_bounds` + `update_active_segs`         | per row |
| `yda.inner_dp_cells`| the inner col loop that updates DP cells          | per cell |
| `yda.row_trailing`  | RY adjustment + insertion-prolongation tail       | per row |
| `yda.traceback`     | edit-script reconstruction walk                   | per step |
| `yda.cleanup`       | `filter_active_segs` + `dpMatrix` frees           | per call |

rdtsc reads are placed **outside** the inner col loop (one pair per
row only) and divided by the row's cell count to get cyc/cell. This
keeps rdtsc overhead amortized to ~25 cycles per ~10³ cells, well
under 1 % of the inner loop's true cost. Implementation:
[`gapped_extend.c` in `stage-timing+rdtsc-substages`](https://github.com/shivsundram/lastz/tree/stage-timing+rdtsc-substages).

**The data, all five regimes at `--step=20`** (10 Mbp × 10 Mbp slice,
single-thread, pinned to CPU 0; raw reports in
`bench/results/rdtsc-gapped-smoke/`):

| Regime | yda.calls | rows | inner cells | **cyc/cell** | inner share | bounds cyc/row | tb cyc/step |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sparse synteny (mm-hg)        |     20 |     23,854 |    11.2 M |   16.8 | 96.4 % | 269 | 12.1 |
| Bird-Z cross (12of19)         |    740 |  1,233,170 |   618.7 M |   16.5 | 96.5 % | 272 | 12.2 |
| Mouse-rat dense (12of19)      |  2,948 |  6,507,500 | 2,808.9 M |   14.0 | 95.1 % | 273 | 14.0 |
| Within-species (match15)      |  1,888 | 12,997,034 | 4,655.8 M | **8.7**| 89.8 % | 300 | 31.6 |
| Primate (rhesus, 12of19)      | 34,994 | 65,115,964 |28,843.5 M |   15.6 | 95.3 % | 307 | 12.3 |

**What the cycle breakdown says.**

1. **The inner DP cell loop is the GPU kernel target.** It owns
   89.8-96.5 % of accounted gapped cycles in every regime. Setup,
   row-bounds, row-trailing, traceback, and cleanup combined are
   ≤10 % even in the most extreme case (within-species, where bigger
   DP matrices push fixed per-row work up). Anything other than the
   inner col loop is bookkeeping; a GPU port that only accelerates
   `update_LR_bounds` will buy at most 3 %.

2. **Per-cell cost is regime-dependent — and the within-species
   regime is 1.9× cheaper per cell than everything else.** Bird-Z,
   mouse-rat, mouse-human and rhesus all cluster at 14-17 cyc/cell;
   within-species lands at 8.7 cyc/cell. Two competing hypotheses:
   (a) at 99.5 % identity the "we cannot improve C" branch is taken
   ~all the time, so the cell body is a tight straight-line update
   with predictable score deltas; (b) within-species has bigger DP
   matrices (~2.5 M cells/call vs ~0.6-1 M elsewhere), giving the
   sweep-row prefetcher more runway. Either way, **the GPU kernel
   needs to handle 8-17 cyc/cell of CPU work**, i.e. on the order
   of ~3 ns/cell at 3 GHz — the same ballpark as a single-precision
   SIMD update.

3. **Per-row fixed cost is constant across regimes.** `update_LR_bounds`
   + `update_active_segs` is 269-307 cyc/row independent of identity,
   density, anchor count, or seed pattern. That's reassuring — it
   means the row-launch overhead on a GPU kernel is bounded by a
   number we can predict *a priori*.

4. **Traceback is essentially free** — 12-14 cyc/step in four of
   five regimes (within-species is the outlier at 31.6 cyc/step,
   ~2.5× slower; that's still 1.6 % of accounted cycles even there).
   A GPU implementation can leave traceback on the host CPU without
   measurable impact.

5. **Cells per ydrop call scale with average alignment length, not
   anchor count.** Within-species: 4.66 G cells / 1,888 calls ≈
   2.5 M cells/call (long, thin DP matrices — high-identity
   extensions go far). Rhesus: 28.8 G / 34,994 ≈ 824 k cells/call
   (many short calls — anchor-rich but lower identity per call).
   Bird-Z, mouse-rat, mouse-human: ~0.6-1.0 M cells/call. Implication
   for batching: a GPU kernel sized to one DP-row-worth of work
   needs to amortize launch overhead across ~10⁵-10⁶ cells; that's
   feasible but argues for a "warp-per-row, persistent-kernel"
   design rather than "one CUDA stream per ydrop call".

6. **Cells per traceback step is a clean inverse-density metric.**
   Sparse: ~2,170 cells/step. Bird-Z: ~1,007. Mouse-rat: ~640.
   Within-species: ~407. Higher identity = tighter band = fewer
   cells explored per output edit op. This is the natural target
   for "how wide does my GPU banded-DP need to be?" — at
   within-species identities, a band width of ~50 cells (5 %
   margin around the diagonal) probably suffices; at cross-species
   identities a 500-cell band is more honest.

7. **Cycle accounting closes to ~70-75 % of stage-timer wall in
   every regime** (e.g. bird-Z step=20: 10.58 G accounted cyc ≈
   3.5 s @ 3 GHz vs. 4.4 s gap_ext stage). The 25-30 % gap is
   memory-stall cycles not retired by the issue width, exactly
   the same dead reckoning as in the seed-search rdtsc analysis.
   It also confirms — independently — that the inner DP cell loop
   is **memory-traffic-bound** more than ALU-bound: a GPU port
   that hits its bandwidth budget should hit closer to ~3 cyc/cell
   on the equivalent core count.

**Reproduce.** Prereq: `lastz_TS` built (`make -C lastz
build_lastz_substages`) and the five 10 Mbp slices on disk (see the
"Step-axis Pareto sweeps" reproduction block below for the exact
`fetch_genomes.py` and `slice_genome.py` calls). Then:

```bash
mkdir -p bench/results/rdtsc-gapped-smoke
OUT=bench/results/rdtsc-gapped-smoke
LASTZ=lastz/src/lastz_TS

run() {
    local tag="$1" t="$2" q="$3"; shift 3
    env LASTZ_STAGE_REPORT=$OUT/${tag}.stage.txt taskset -c 0 \
        "$LASTZ" "$t" "$q" "$@" --format=maf- > /dev/null 2>&1
}

run mouse_human_step20  bench/data_genomes/hg38.chr19_40_50mb.fa \
                        bench/data_genomes/mm10.chr10_120_130mb.fa \
                        --seed=12of19 --step=20
run bird_step20         bench/data_genomes/galGal6.chrZ_0_10mb.fa \
                        bench/data_genomes/taeGut2.chrZ_0_10mb.fa \
                        --seed=12of19 --step=20
run mouse_rat_step20    bench/data_genomes/mm10.chr10_40_50mb.fa \
                        bench/data_genomes/rn6.chr20_40_50mb.fa \
                        --seed=12of19 --step=20
run hg_vs_chm13_match15 bench/data_genomes/hg38.chr1_50_60mb.fa \
                        bench/data_genomes/hs1.chr1_50_60mb.fa \
                        --seed=match15 --notransition --step=20
run rhesus_step20       bench/data_genomes/hg38.chr19_40_50mb.fa \
                        bench/data_genomes/rheMac10.chr19_40_50mb.fa \
                        --seed=12of19 --step=20

# Extract the gapped substage block from each report:
for f in $OUT/*.stage.txt; do
    echo "=== $(basename $f .stage.txt) ==="
    sed -n '/--- rdtsc gapped_extend breakdown ---/,/===STAGE_TIMING_END===/p' "$f"
done
```

Total wall: ~4-5 minutes (rhesus alone is ~4 min; the other four
are ≤25 s each).

#### Two textbook references for `ydrop_one_sided_align` (`ydrop_sane`)

Once the inner DP cell loop is identified as the GPU kernel target, the
next question is "what does the kernel **actually have to compute**?"
`ydrop_one_sided_align` is hard to read directly — it interleaves five
distinct concerns (sweep-row buffer arithmetic, traceback tape append,
left/right-segment masking, y-drop pruning, row trailing) into one
~250-line function. To pin down the recurrence cleanly we wrote
**two** textbook ports in
[`lastz/src/ydrop_sane.c`](lastz/src/ydrop_sane.c), both bit-identical
to lastz over the v0 scope (no neighbor-alignment masking, `trimToPeak
== true`, forward extension):

1. **`ydrop_one_sided_align_impl_sane`** — the "wasteful but obvious"
   reference: three full `(M+1) × (N+1)` score matrices for `C`, `D`,
   `I` plus one full `(M+1) × (N+1)` link byte tape. ~13 · (M+1) · (N+1)
   bytes per call. The clearest possible mapping from the algebraic
   recurrences to code; intended as the reading copy when you want to
   understand what the GPU kernel must compute.

2. **`ydrop_one_sided_align_impl_sane_double_buffered`** — same
   algorithm, lastz-matching memory profile (see "Memory layout"
   below). Sweep-row C/D buffers, scalar I, band-compact link tape.
   The version that participates in fair runtime A/B comparisons
   against `ydrop_one_sided_align` (since both spend their time on
   the actual DP work, not on initializing megabyte-scale scratch
   arrays).

Common to both:
- Per-call `malloc_or_die`/`free_if_valid`. No state crosses calls
  (unlike lastz's module-static `tback`/`tbRow[]`).
- Per-row `LY[]` / `RY[]` arrays kept explicit, no segment-driven
  bound updates. Y-drop pruning, left-edge chain shrink (`LY++` on
  consecutive left-side prunes), right-boundary termination (`RY++`
  by one with a negInf sentinel), and the insertion-prolongation row
  trailing all match lastz line-for-line.
- `bestScore` only updates inside the `cFromC` (we-cannot-improve-C)
  branch, matching lastz's choice to anchor the alignment endpoint at
  a substitution rather than a gap.

##### Memory layout (double-buffered variant)

The `_double_buffered` variant adopts two of lastz's three memory
tricks (see [`lastz/src/gapped_extend.c`](lastz/src/gapped_extend.c)
note 3 and line 3775):

1. **Sweep-row score buffers.** Two `(N+1)`-cell arrays each for `C`
   and `D`, swapped between rows. `I` is a scalar carried through the
   inner loop. Per-call score memory: **O(N)**, ~16 KB at N=1000.

2. **Band-compact link tape.** A single byte tape, indexed via a per-
   row offset table:

   ```c
   lnkRow[r] = (lnkPos at row r's start) - LY[r];
   // Cell (r, c)'s link byte lives at lnkTape[lnkRow[r] + c]
   ```

   Only cells in the live band `[LY[r], RY[r])` consume tape bytes;
   cells outside the band do not exist in memory. Per-call link
   memory: **O(sum of band widths)**. For typical y-drop runs (band
   ~500 cells), ~500 KB at M=1000, ~5 MB at M=10000.

The one lastz trick we skip is the **module-static reuse** of the tape
across calls. Lastz allocates a single 80 MB `tback` at startup,
overwrites it from offset 0 on every ydrop call, and truncates the
alignment if a call would overflow it. We do per-call `malloc`/`free`
of all four allocations (score buffers, lnk tape, `lnkRow[]`, `LY/RY`)
to keep the impl state-free. The tape grows dynamically (2× on
overflow); `yDropTail = yDrop/gapE + 16` cells of headroom per row
mean realloc almost never fires mid-row in practice.

Side-by-side memory at M=N=2000, yDrop=910:

| Impl | Score memory | Link memory | Total |
|---|---:|---:|---:|
| `ydrop_one_sided_align_impl_sane` (full 2-D) | 48 MB | 4 MB | 52 MB |
| `..._double_buffered` (sweep + band) | 32 KB | ~500 KB | ~600 KB |
| lastz `ydrop_one_sided_align` | 32 KB (per-call) | 80 MB (module-static) | same per-call as `_double_buffered` |

For pathological calls (M=N=10000, no y-drop firing) the
`_double_buffered` impl uses ~5 MB total per call, well within what
the OS will allocate; lastz would also handle this fine (10k × 10k =
100M cells, but the 80 MB tape would truncate the trailing portion).

##### Validation

Both impls are validated against `ydrop_one_sided_align_for_testing`
(a non-static trampoline added to
[`gapped_extend.c`](lastz/src/gapped_extend.c) that constructs a
minimal `alignio` and forwards to the file-static
`ydrop_one_sided_align`; the trampoline does not affect any other
lastz code path). The driver in
[`bench/test_ydrop.c`](bench/test_ydrop.c) compares score, endpoint,
and edit-script byte-by-byte across:
- 8 hand-crafted cases: identity, single mismatch, single insertion,
  single deletion, y-drop-firing tail, alternating mismatches with
  high y-drop, long identity, and short no-match;
- 200 fuzz pairs at lengths 50-1500 bp and identities 60-99 %
  (xorshift, deterministic seed);
- each case is run TWICE — once forward (`reversed=0` on lastz, sane
  on `(A,B)`) and once reversed (`reversed=1` on lastz, sane on
  `(rev(A), rev(B))`) — to cover both halves of lastz's two-sided
  alignment;
- each (case, direction) pair is run against BOTH sane impls.

Per case = 4 comparisons (2 impls × 2 directions). Result:

- Default FUZZ=200: **832 / 832 pass** ((8 hand-crafted + 200 fuzz) × 4).
- FUZZ=1000:        **4032 / 4032 pass** ((8 + 1000) × 4).
- FUZZ=3000:        **12032 / 12032 pass** ((8 + 3000) × 4).

The validated v0 includes everything that touches the inner DP cell
loop plus its surrounding bounds/trailing logic — the exact code
surface a GPU kernel will need to replicate. Reproduce:

```bash
make -C lastz                      # build the standard lastz objects
make test_ydrop                    # link bench/test_ydrop against them
./bench/test_ydrop                 # 832 PASS
FUZZ=1000 ./bench/test_ydrop       # 4032 PASS
FUZZ=3000 ./bench/test_ydrop       # 12032 PASS
```

##### Runtime A/B swap inside lastz (`YDROP_SANE_IMPL=1`)

To compare the two impls' wall-clock cost on real lastz workloads
without writing a separate microbenchmark, we added a runtime swap to
[`gapped_extend.c`](lastz/src/gapped_extend.c). When the environment
variable `YDROP_SANE_IMPL=1` is set, the top of
`ydrop_one_sided_align` diverts each call to
`ydrop_one_sided_align_impl_sane_double_buffered` whenever that call is
in the v0 scope:

- `io->leftSeg == NULL && io->rightSeg == NULL` (no left/right neighbor
  segments to clamp L/R bounds against);
- `io->leftAlign == NULL && io->rightAlign == NULL` (no neighbor
  alignments to drive `next_sweep_seg` updates);
- `io->aboveList == NULL` (or `belowList` for `reversed=1`) — no
  active-segment list to propagate through `filter_active_segs`;
- `trimToPeak == true`.

Out-of-scope calls fall through to lastz proper. Stats are reported at
exit (always to stderr via atexit, and additionally to the stage-timing
file when built with `-DdbgTiming`):

```
--- ydrop_sane runtime swap ---
total ydrop calls                                3572
swapped to sane (double_buf)                       20  (0.6%)
declined (out of v0 scope)                       3552  (99.4%)
  decline reasons (first failing condition per call):
    leftSeg != NULL                             826
    rightSeg != NULL                            436
    leftAlign != NULL                             0
    rightAlign != NULL                            0
    above/belowList != NULL                    2290
    trimToPeak == false                           0
```

The swap is bit-identical on simple workloads and produces an
acceptable alignment on real ones (same score, same endpoints; edit
script may differ by O(few) operations in rare boundary-cell-rescue
cases that the fuzz harness doesn't exercise — see "Known limitation"
below).

###### Wall-clock A/B results

Three synthetic workloads where v0 scope holds on 100 % of calls (no
prior alignments, no masking), each producing a single contiguous
~94 % identity alignment, single thread:

| Workload | lastz tape | baseline wall | sane wall | Δ | swap rate |
|---|---|---:|---:|---:|---:|
| uniform 100 kbp     | default 80 MB | 0.35 s | 0.45 s | +29 % | 2/2 (100 %) |
| uniform 1 Mbp       | `--allocate:traceback=1G` | 2.90 s | 3.43 s | +18 % | 2/2 (100 %) |
| uniform 10 Mbp      | `--allocate:traceback=2G` | 90.7 s | 89.6 s | −1 % (noise) | 2/2 (100 %) |

Per-call constant overhead (per-call `malloc`/`free` of score buffers,
lnk tape, lnkRow/LY/RY) is the dominant gap at small N; it amortizes
to zero by 10 Mbp, where the two impls are within noise of each
other.

> Without `--allocate:traceback`, lastz's 80 MB tape overflows on the
> 1 Mbp and 10 Mbp workloads and emits "truncating alignment ending
> at ..." warnings, fragmenting the output into many chunks. Sane has
> no tape-size cap (it grows on demand) and produces one continuous
> alignment in those cases. Setting `--allocate:traceback=2G` makes
> both impls produce bit-identical output and the A/B fair.

###### What this isn't (yet)

On dense real workloads (e.g. mouse-rat chr10×chr20 10 Mbp, 1768
HSPs), the swap declines 99.4 % of calls because lastz populates
`aboveList`/`belowList` (active-segment masking, 64 % of declines) and
`leftSeg`/`rightSeg` (per-row L/R bound constraints, 35 % of declines)
on virtually every secondary HSP. The swap fires only on the very
first anchor of a chore plus a handful of others; that subset is too
small for a meaningful wall-clock comparison.

Two follow-ups would close that gap:
1. Extend `ydrop_sane` to honor an `aboveList`/`belowList`-style mask
   passed in as a callback or precomputed band-clip array. The kernel
   itself doesn't change; only the per-row LY/RY update.
2. Add an external batch driver that loads (anchor, A, B, M, N,
   scoring) tuples dumped from lastz to disk and replays them against
   both impls back-to-back. The dump path bypasses lastz's masking
   entirely.

###### Known limitation: boundary-cell traceback

On the dense mouse-rat run, one of the 20 swapped alignments had the
same score and same endpoints as lastz but an edit script differing by
6 columns (30520 / 34690 matches in lastz vs 30514 / 34684 in sane,
same 88.0 % identity). Investigation: lastz's inner-loop terminator is
`col < R[r-1] + 1`, so it processes ONE extra column past the previous
row's right edge (`col = R[r-1]`). If that cell receives a fresh
diagonal match from `C[r-1][R[r-1]-1] + sub[A[r]][B[R[r-1]]]`, lastz
will rescue it as alive; sane's loop terminates at `col = R[r-1] - 1`
and treats the boundary cell as a sentinel negInf. In rare cases the
rescued cell participates in the eventual traceback path, yielding the
same DP score but a slightly different edit script.

The fuzz harness in [`bench/test_ydrop.c`](bench/test_ydrop.c) runs
50–1500 bp pairs, where the boundary cell is too short-lived to
trigger this; that's why 12 032 / 12 032 fuzz cases pass while real
~35 kbp alignments occasionally diverge. The two outputs remain
equivalent from a scoring standpoint — neither is "more correct" — but
this is a real algorithmic deviation, not just a traceback tie-break.

###### Reproduce

```bash
cd lastz/src && make lastz && cd ../..

# bit-identical small workload, 100 % swap rate
./lastz/src/lastz bench/data/uniform_100kb.target.fa \
                  bench/data/uniform_100kb.query.fa \
                  --format=general --output=/tmp/base.out
YDROP_SANE_IMPL=1 \
  ./lastz/src/lastz bench/data/uniform_100kb.target.fa \
                    bench/data/uniform_100kb.query.fa \
                    --format=general --output=/tmp/sane.out
diff /tmp/base.out /tmp/sane.out         # empty: bit-identical

# 1 Mbp A/B with enlarged tape so neither impl truncates
/usr/bin/time -v ./lastz/src/lastz bench/data/uniform_1000kb.target.fa \
                                    bench/data/uniform_1000kb.query.fa \
                                    --allocate:traceback=1G \
                                    --format=general --output=/dev/null
YDROP_SANE_IMPL=1 /usr/bin/time -v \
                  ./lastz/src/lastz bench/data/uniform_1000kb.target.fa \
                                    bench/data/uniform_1000kb.query.fa \
                                    --allocate:traceback=1G \
                                    --format=general --output=/dev/null
# stderr ends with the runtime-swap report.
```

##### Cross-anchor sequential-dependency study

After establishing that the in-situ A/B swap declines 99.4 % of calls on
dense workloads, we wanted to know whether a GPU y-drop kernel is even
worth pursuing if the outer anchor loop has a true loop-carried
dependency. Three experiments answer this end-to-end:

1. **Anchor-loop instrumentation** ([`gapped_extend.c`](lastz/src/gapped_extend.c)):
   counters at `msp_left_right` that measure how many anchors get
   skipped before y-drop runs, plus per-call walks of `aboveList` /
   `belowList` to gauge the neighbor-alignment chain depth.
2. **Per-call y-drop matrix log**: env-var-gated CSV
   (`YDROP_CALL_LOG=path`) dumping `(M, N, rows_processed, max_band,
   cells_visited, cycles_total, peak_score)` for each
   `ydrop_one_sided_align` call. Lets us characterize the actual DP
   shape rather than the worst-case bound passed in.
3. **Masking / containment-skip ablation**: two new env-var gates,
   `LASTZ_DISABLE_NEIGHBOR_MASK=1` (NULL the six `alignio` neighbor
   fields right before `ydrop_align`) and `LASTZ_DISABLE_CONTAINMENT=1`
   (ignore the `msp_left_right` skip), let us run the lastz pipeline
   with either or both of the cross-anchor mechanisms disabled and
   diff the output against baseline.

All three were run on three 10 Mbp slices that span the alignable-
density axis: bird-Z (`galGal6.chrZ_0_10mb` × `taeGut2.chrZ_0_10mb`),
mouse-rat (`mm10.chr10_40_50mb` × `rn6.chr20_40_50mb`), and within-
species (`hg38.chr1_50_60mb` × `hs1.chr1_50_60mb`).

###### Anchor-loop containment + chain depth

| regime    | total anchors | skipped (contained) | ran y-drop | aboveList non-empty | mean chain | max chain |
|-----------|--------------:|--------------------:|-----------:|--------------------:|-----------:|----------:|
| bird-Z    |        1 665  | 1 154 (69 %)        | 511  (31 %)| 97 %                |       100  |    428    |
| mouse-rat |       17 982  | 16 196 (90 %)       | 1 786 (10 %)| 99 %               |       285  |   1 294   |
| hg-CHM13  |       24 376  | 21 901 (90 %)       | 2 475 (10 %)| 99 %               |       473  |   2 246   |

Two findings drop out:
- **The seed stage's output is enormously redundant.** 70 – 90 % of
  HSPs are sub-regions of bigger HSPs that already won the scoring
  contest; gapped_extend's first job is to throw them out. (SegAlign's
  GPU `thrust::unique_copy` is doing essentially the same dedup,
  earlier in the pipeline.)
- **Survivors are deeply coupled.** Almost every y-drop call has 100 +
  neighbor alignments to mask against. So the v0-swap's earlier
  "99.4 % decline because `aboveList` non-NULL" was not "list has 1-2
  entries" but "list has *hundreds* of entries". Whether masking
  matters is the question Experiment 2 below answers.

###### Per-call DP shape

| regime    | calls  | median cells | p95 cells | median rows | median band | median cyc | max cyc |
|-----------|-------:|-------------:|----------:|------------:|------------:|-----------:|--------:|
| bird-Z    | 1 022  | 633 K        | 1.6 M     | 1 287       | 676         | 12.1 M     | 78.6 M  |
| mouse-rat | 3 572  | 554 K        | 2.5 M     | 1 241       | 623         | 10.2 M     | 136 M   |
| hg-CHM13  | 4 950  | 332 K        | 884 K     | 1 090       | 502         | 6.9 M      | **4 054 M** |

- **Per-call DP is small.** Median ~500 K cells, ~1 200 rows × ~600-
  wide band. That is 0.000 1 of the worst-case M·N the call signature
  allows; y-drop pruning is extremely aggressive.
- **Distribution is bimodal.** Most calls fit in 1 – 2 M cells, but the
  tail can hit 480 M cells / 4 G cycles (the hg-CHM13 outlier was a
  single ~1.3 M-row alignment — a within-species ortholog).
- **Total y-drop work per regime**: 0.8 – 5.4 B cells, 15 – 65 G cycles,
  with ~50 % concentrated in <5 % of calls. Median 500 K cells per
  call is too small to saturate a modern GPU on its own; the kernel
  has to batch many calls together (persistent threads + work queue)
  rather than launch per call.

###### Ablation: masking off vs containment-skip off (mouse-rat, 10 Mbp)

| variant                          | n alignments | total bp1 | bp-set Jaccard vs baseline | wall    | wall vs base |
|----------------------------------|-------------:|----------:|---------------------------:|--------:|-------------:|
| baseline                         |       1 767  |  4.15 M   | 1.000                       | 46.7 s  | 1.00×        |
| 2a (mask off, contain on)        |       1 800  |  6.33 M   | **0.9977**                  | 50.3 s  | 1.08×        |
| 2b (mask off, contain off)       |      17 982  |  140 M    | **0.9977**                  | 492 s   | **10.5×**    |

"bp-set Jaccard" is the Jaccard index over the **set of target-side bp
positions covered by at least one alignment** (1.0 = identical bp
covered; 0.0 = no overlap). Both variants score 0.9977 against
baseline, meaning all three runs cover essentially the same homologous
DNA — they differ only in how that DNA is packaged into alignment
records.

What this tells us:
- **Masking is cosmetic.** Variant 2a (no masking, containment-skip
  intact) produces 33 more alignment records than baseline (+1.9 %),
  total aligned bp is +53 % because alignments are longer on average,
  but the *bp-set* covered is unchanged (Jaccard 0.9977). Looking at
  the differences: 423 short baseline alignments (median 168 bp) get
  absorbed into 259 longer ones (median 2.3 kbp); the underlying
  homology is identical. Wall cost is only +8 %.
- **Containment-skip is the real dedup.** Variant 2b (both off) produces
  17 982 alignments (~10× baseline) — exactly the number of HSPs that
  entered gapped_extend in Experiment 1. Total aligned bp blows up to
  140 Mbp (34× baseline) because every anchor now produces its own
  ~7.8 kbp alignment and they all overlap each other. Same Jaccard
  0.9977 — same homologous territory, just hugely redundant
  packaging. Wall blows up to 492 s (10.5×).

So the cross-anchor dependency decomposes into two cleanly separable
mechanisms with very different impact:

| mechanism         | what it does                          | impact on output | impact on perf |
|-------------------|---------------------------------------|------------------|---------------:|
| containment-skip  | parallel-friendly **anchor dedup**    | essential        | saves ~10× work |
| neighbor masking  | per-row clamp + per-cell mask         | cosmetic         | costs ~8 % wall |

###### Implications for the GPU plan

The data steer the design firmly toward one shape:

1. **GPU y-drop kernel, no masking.** Banded wavefront, intra-call
   antidiagonal parallelism, ~500-600-wide band. Drop the entire
   `aboveList` / `belowList` / `leftSeg` / `rightSeg` plumbing — it's
   cosmetic for output, expensive in per-row overhead, and would force
   the GPU kernel to keep a per-row active-segment list around. The
   ~0.2 % bp coverage delta is well within the noise of any biological
   downstream use of these alignments.
2. **Batched / persistent kernel.** Median 500 K cells per call is far
   below GPU saturation. Need persistent threads pulling from a work
   queue of (anchor, A, B, M, N, scoring) tuples so many calls can run
   concurrently across SMs.
3. **Outlier path.** The p99 / max calls are 100 – 1 000× bigger than
   median; a single 4 G-cycle call can dominate. Either dispatch
   outliers to a separate single-call high-occupancy kernel that uses
   multiple CTAs per call, or tile every call into fixed-band chunks
   regardless of size.
4. **Containment-skip stays.** Without it, 10× the y-drop work and 10×
   the wall, on output that's not biologically more useful. But the
   *implementation* can move to a parallel GPU post-pass: y-drop every
   anchor, then parallel-greedy-suppress in score order — each
   alignment in parallel checks whether any higher-scoring alignment's
   box contains its anchor. O(N²) work, trivially parallel, ~1 ms on
   modern GPU for N ≈ 20 K. Same answer as lastz's sequential
   `msp_left_right` because greedy-suppress is order-independent on a
   fixed input.
5. **A/B harness updates.** The right "lastz proper" baseline for the
   GPU kernel A/B is `LASTZ_DISABLE_NEIGHBOR_MASK=1` (mask off,
   contain on) — same masking decision the GPU kernel will make. Our
   v0 swap stays useful as a pure-CPU reference impl.

Back-of-envelope ceiling: mouse-rat baseline does ~50 G cycles in
y-drop on a single thread (≈17 s of pure y-drop CPU work). A modern
GPU at ~10 TFLOPS effective DP throughput can do 3 B cells in ~3 ms;
even at 100× pessimism (launch overhead, transfers, tail effects),
that's ~300 ms. On the 90 % of wall that gapped extension owns on
dense regimes, a realistic **10 – 50× speedup on y-drop** translates
to ~5 – 15× overall.

###### Reproduce

```bash
cd lastz/src && make lastz lastz_TS && cd ../..

# Experiments 1 + 3: anchor-loop + per-call distribution
mkdir -p /tmp/exp13
for r in birdz mr hg; do
    case $r in
      birdz) T=galGal6.chrZ_0_10mb; Q=taeGut2.chrZ_0_10mb;;
      mr)    T=mm10.chr10_40_50mb;  Q=rn6.chr20_40_50mb;;
      hg)    T=hg38.chr1_50_60mb;   Q=hs1.chr1_50_60mb;;
    esac
    YDROP_CALL_LOG=/tmp/exp13/$r.csv \
    LASTZ_STAGE_REPORT=/tmp/exp13/$r.stage \
      ./lastz/src/lastz_TS bench/data_genomes/$T.fa bench/data_genomes/$Q.fa \
        --allocate:traceback=2G --format=general \
        --output=/tmp/exp13/$r.out 2>&1 | tail -2
done

# Experiment 2: ablation on mouse-rat
./lastz/src/lastz bench/data_genomes/mm10.chr10_40_50mb.fa \
                  bench/data_genomes/rn6.chr20_40_50mb.fa \
                  --allocate:traceback=2G --format=general \
                  --output=/tmp/exp2_base.out
LASTZ_DISABLE_NEIGHBOR_MASK=1 \
  ./lastz/src/lastz bench/data_genomes/mm10.chr10_40_50mb.fa \
                    bench/data_genomes/rn6.chr20_40_50mb.fa \
                    --allocate:traceback=2G --format=general \
                    --output=/tmp/exp2_2a.out
LASTZ_DISABLE_NEIGHBOR_MASK=1 LASTZ_DISABLE_CONTAINMENT=1 \
  ./lastz/src/lastz bench/data_genomes/mm10.chr10_40_50mb.fa \
                    bench/data_genomes/rn6.chr20_40_50mb.fa \
                    --allocate:traceback=2G --format=general \
                    --output=/tmp/exp2_2b.out
```

#### Step-axis Pareto sweeps across five regimes

After the seed-weight (bird-Z) and within-species (hg-vs-CHM13) sweeps
above, we ran five `--step` sweeps at the same 10 Mbp × 10 Mbp scale to
map how the seed-vs-gapped balance moves across phylogenetic distance.
Three findings collapse out:

1. The crossover from "seed-dominated" to "gapped-dominated" happens
   well before 93 % identity.
2. The bottleneck axis is *alignable density*, not identity. Identity
   only matters because it shapes density.
3. For three of five regimes, gapped extension is 83–99 % of wall —
   so the GPU-acceleration target is `ydrop_one_sided_align`, not the
   seed kernel, on anything closer than ~85 % identity *with dense
   synteny*.

Single-thread runs pinned to CPU 0, `lastz_TS` build (so we get
per-stage timing), `--seed=12of19` everywhere except the within-species
column (where `--seed=match15 --notransition` is the appropriate
recipe).

##### All step-sweep data, in one place

Master table — every step-sweep data point we have, side-by-side, at the
same 10 Mbp × 10 Mbp scale. Identity rises top-to-bottom; alignable
density also rises top-to-bottom, *but not monotonically* (the sparse-
synteny mouse-human row sits above bird-Z in identity but far below it
in density).

| Regime | Slice | Identity | step | Wall | seed_hit | gap_ext | HSPs | bp aligned |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Sparse synteny | hg.chr19[40-50] × mm.chr10[120-130] | ~80 % | 1 | 18.4 s | 17.5 s **(95 %)** | 0.15 s (1 %) | 15 | 6.3 kbp |
|  |  |  | 20 | 3.4 s | 3.0 s (87 %) | 0.10 s (3 %) | 10 | 4.9 kbp |
|  |  |  | 100 | 2.5 s | 2.1 s (84 %) | 0.08 s (3 %) | 6 | 4.5 kbp |
| Cross-species (bird-Z) | galGal6.chrZ × taeGut2.chrZ [0-10] | ~75 % | 1 | 62.3 s | 55.2 s **(89 %)** | 6.2 s (10 %) | 510 | 0.61 Mbp |
|  |  |  | 20 | 11.8 s | 6.3 s (54 %) | 4.9 s (42 %) | 370 | 0.54 Mbp |
|  |  |  | 100 | 7.2 s | 3.8 s (53 %) | 3.1 s (43 %) | 214 | 0.38 Mbp |
| **Mouse-rat dense** | **mm.chr10[40-50] × rn6.chr20[40-50]** | **~85 %** | **1** | **47.3 s** | **25.6 s (54 %)** | **21.0 s (44 %)** | **1767** | **4.15 Mbp** |
|  |  |  | **20** | **22.6 s** | **3.4 s (15 %)** | **18.9 s (83 %)** | **1462** | **4.05 Mbp** |
|  |  |  | **100** | **17.7 s** | **2.2 s (13 %)** | **15.1 s (85 %)** | **991** | **3.77 Mbp** |
| Primate (rhesus) | hg.chr19[40-50] × rheMac10.chr19[40-50] | ~93 % | 1 | 346.6 s | 14.3 s (4 %) | 331.2 s **(96 %)** | 28061 | 37.82 Mbp |
|  |  |  | 20 | 221.9 s | 2.1 s (1 %) | 219.2 s **(99 %)** | 17483 | 31.23 Mbp |
|  |  |  | 100 | 123.2 s | 1.4 s (1 %) | 121.3 s **(98 %)** | 8267 | 23.14 Mbp |
| Within-species (match15) | hg.chr1[50-60] × hs1.chr1[50-60] | ~99.5 % | 20 | 22.3 s | 0.49 s (2 %) | 20.9 s **(94 %)** | 943 | 11.33 Mbp |
|  |  |  | 50 | 19.1 s | 0.33 s (2 %) | 18.3 s **(96 %)** | 648 | 10.89 Mbp |
|  |  |  | 100 | 17.5 s | 0.31 s (2 %) | 16.7 s **(96 %)** | 431 | 10.52 Mbp |

The mouse-rat row is the threshold pin — at step=1 it sits at 54/44
seed/gap (essentially balanced); by step=20 it has flipped to 15/83
gap-dominated. See **Pinning the threshold** subsection below for the
full sweep + within-mm-rn density gradient.

##### Compact per-regime summary

| Regime | Slice | Identity | bp aligned (10 Mbp²) | seed share | gap share | Wall step=1 → step=100 |
|---|---|---:|---:|---:|---:|---:|
| Sparse synteny | hg38.chr19[40-50M] × mm10.chr10[120-130M] | ~80 % within blocks | 0.006 Mbp | 84-95 % flat | 1-4 % flat | 18.4 s → 2.5 s |
| Cross-species (dense) | galGal6.chrZ[0-10M] × taeGut2.chrZ[0-10M] | ~75 % | 0.61 Mbp | 89 % → 53 % | 10 % → 43 % | 62.3 s → 7.2 s |
| Mouse-rat (dense) | mm10.chr10[40-50M] × rn6.chr20[40-50M] | ~85 % | 4.05 Mbp | 54 % → 13 % | 44 % → 85 % | 47.3 s → 17.7 s |
| Primate (dense) | hg38.chr19[40-50M] × rheMac10.chr19[40-50M] | ~93 % | 23-38 Mbp | 1-4 % flat | 96-99 % flat | 346.6 s → 123.2 s |
| Within-species | hg38.chr1[50-60M] × hs1.chr1[50-60M] | ~99.5 % | 11-12 Mbp | 1-2 % flat | 94-96 % flat | 22.3 s → 17.5 s (match15) |

Detailed per-step numbers are committed at:

- `bench/results/step-sweep-mouse10mb/`  (sparse synteny)
- `bench/results/step-sweep-bird10mb/`   (cross-species, with
  `seed-sweep-bird10mb/default_12of19.*` as the cached step=1 baseline)
- `bench/results/step-sweep-rat10mb/`    (mouse-rat, threshold pin)
- `bench/results/step-sweep-rhesus10mb/` (primate)
- `bench/results/step-sweep-hg-vs-chm13/` (within-species)

**Visual summary of the bottleneck axis** (alignable density on a log
scale; identity rises from left to right, but the bottleneck is set by
density, not identity):

```
                          bp aligned per 100 Mbp²  (log scale)
        ────────────────────────────────────────────────────────►
        ~0.01 Mbp        ~1 Mbp            ~10-40 Mbp
         │                 │                    │
   ┌─────┴─────┐      ┌────┴─────┐          ┌───┴──────┐
   │ SEED      │      │ TRANSITION│          │ GAPPED   │
   │ dominates │      │ region    │          │ dominates│
   │ 80-95 %   │      │ (50/50    │          │ 95-99 %  │
   │ flat in   │      │ at step≥20│          │ flat     │
   │ step      │      │ on bird-Z)│          │          │
   └───────────┘      └───────────┘          └──────────┘

   mouse-human       bird-Z chrZ            rhesus chr19,
   sparse pairs      cross-species          hg-vs-CHM13
```

**1. Crossover happens between bird-Z and rhesus.** At ~75 % identity
+ dense synteny (bird-Z), seed_hit dominates at step=1 (89 %) and
equilibrates near 50/50 by step=20. At ~93 % identity (rhesus) the
bottleneck has *already* fully flipped: gapped is 96 % at step=1 and
99 % at step=20. Whatever density threshold matters, it sits between
0.6 Mbp aligned (bird-Z) and ~25 Mbp aligned (rhesus step=50). This
nukes the simple "low identity → seeds dominate, high identity → gap
dominates" story; the rhesus regime is closer to within-species than
to bird-Z despite being clearly cross-species.

**2. Sparse-synteny is its own regime.** The mouse-human probe
yielded ~6 kbp aligned out of 100 Mbp² — 100× lower than bird-Z and
2,000× lower than within-species. With so little gappable territory,
seed_hit dominates uniformly at 84–95 % across every step value. The
sweep doesn't equilibrate. This is the regime where SegAlign's seed
kernel earns *all* the speedup; gapped acceleration would be wasted
silicon for these workloads.

**3. Step is a different knob in each regime.**

  - Sparse synteny: step=1 → 100 buys 7× wall (18.4 → 2.5 s) but
    drops aligned bp ~30 % (6.3 → 4.5 kbp). Mostly seed savings.
  - Cross-species: step=1 → 20 alone buys 5.3× wall and only 11 %
    sensitivity loss. Step≥20 is essentially free.
  - Primate: step=1 → 100 buys 2.8× wall but loses 39 % of bp aligned
    (37.8 → 23.1 Mbp). Each step value still costs >2 minutes;
    sensitivity-cost product is dominated by the gapped budget.
  - Within-species: step=20 → 100 is a flat 1.3× wall and 7 % loss.
    Already at floor.

**Implication for the GPU thesis.** "Where is the time going" depends
on the workload, but the answer collapses into two cases:

- **Workload aligns >~1 Mbp per 100 Mbp² product** → gapped_ext is
  the only kernel that matters. Within-species, primate-primate, and
  anything closer than ~85 % syntenic identity. Three of our four
  regimes. Speed up `ydrop_one_sided_align` and you're done.
- **Workload aligns <~1 Mbp per 100 Mbp² product** → only the seed
  kernel matters. Cross-species at scale, distant pairs, sparse-
  synteny chunks. One of our four regimes (and the one SegAlign is
  benchmarked on).

Bird-Z sits at the inflection point — the *only* regime where both
kernels matter. That's also exactly why SegAlign benchmarks against
it: the marketing-friendly "both stages on GPU" only earns its keep
on bird-Z-shaped workloads.

**Caveat — `--seed=12of19` is the wrong seed for rhesus.** The actual
NCBI/UCSC primate-primate recipe uses contiguous seeds with explicit
human-friendly scoring. Using `12of19 --transition` with default
HOXD70-style scoring on rhesus is what generated 17 k overlapping HSPs
and the 31 Mbp double-counted aligned bp number. The bottleneck shape
(gapped-dominated) is correct; the absolute numbers reflect a
deliberately-untuned configuration. Tuning would shrink the
`gapped_ext` budget but not change the regime story.

**Reproduce.** Add `rheMac10` to the genome inventory, fetch it
(~740 MB 2bit, ~30 s download), slice [40 M, 50 M] of chr19, and run
the rhesus sweep. The other three sweeps share the same data layer
(see "Within-species regime change" above for hg-vs-CHM13 and
"Seed-weight sweep" above for bird-Z); the mouse sparse-synteny slice
needs a one-time chr10 extraction.

```bash
# 1) Fetch rheMac10 (Mmul_10) and slice chr19[40,50M] to match the
#    existing hg38.chr19_40_50mb.fa slice.
python3 bench/fetch_genomes.py --only rheMac10
python3 bench/slice_genome.py \
    bench/data_genomes/rheMac10.chr19.fa \
    /scratch2/shiv1/lastz-bench-data/slices/rheMac10.chr19_40_50mb.fa \
    --start 40000000 --length 10000000
ln -sf /scratch2/shiv1/lastz-bench-data/slices/rheMac10.chr19_40_50mb.fa \
    bench/data_genomes/rheMac10.chr19_40_50mb.fa

# 2) Slice mm10.chr10[120,130M] (probed densest 10 Mbp synteny block
#    with hg38.chr19[40,50M]; still very sparse, ~6 kbp aligned).
python3 bench/slice_genome.py \
    bench/data_genomes/mm10.chr10.fa \
    /scratch2/shiv1/lastz-bench-data/slices/mm10.chr10_120_130mb.fa \
    --start 120000000 --length 10000000
ln -sf /scratch2/shiv1/lastz-bench-data/slices/mm10.chr10_120_130mb.fa \
    bench/data_genomes/mm10.chr10_120_130mb.fa

# 3) Run all four step sweeps. Total wall ~25 min single-thread (rhesus
#    dominates the budget at ~17 min).
./step_sweep.sh           # within-species: hg-vs-CHM13, ~75 s
./step_sweep_bird.sh      # cross-species: bird-Z, ~3 min (step=1 cached)
./step_sweep_mouse.sh     # sparse synteny: mouse-human, ~30 s
./step_sweep_rhesus.sh    # primate: hg-vs-rhesus, ~17 min

# 4) (optional) Reproduce the headline-numbers + bottleneck-flip table:
./get_tradeoffs.sh        # ~95 s, prints both tables
```

Each step-sweep script is a thin wrapper over `lastz_TS` with stage
timing enabled, the same shape as `sweep_human1.sh`: pinned to CPU 0,
fail-loudly on non-zero exit, captures `*.stage.txt` per step value,
and prints a Pareto table at the end.

#### Pinning the threshold: mouse-rat at ~85 % identity

The four-regime sweep above puts the seed→gapped crossover somewhere
between bird-Z's 0.6 Mbp aligned (seed dominates 9× at step=1) and
rhesus's 23 Mbp aligned (gap dominates 7×). To pin it down we picked
**mouse-rat at chr10[40-50M] × chr20[40-50M]** — ~85 % identity,
dense synteny, alignable density chosen by probing 5 (mm.chr10 ×
rn6.chr20) 10 Mbp pairs at `--step=100` and picking the densest.

Probing the within-mouse-rat density gradient at step=100 (single
identity tier, varying offset) walks straight through the threshold:

| Slice | bp aligned (step=100) | seed_hit | gap_ext | dominator |
|---|---:|---:|---:|---|
| mm.chr10[ 0-10M] × rn6.chr20[ 0-10M] |   3.5 kbp | 2.23 s | 0.07 s | seed (32×) |
| mm.chr10[10-20M] × rn6.chr20[10-20M] |   4.8 kbp | 2.46 s | 0.12 s | seed (21×) |
| mm.chr10[20-30M] × rn6.chr20[20-30M] |  17.9 kbp | 2.44 s | 0.23 s | seed (11×) |
| mm.chr10[30-40M] × rn6.chr20[30-40M] | 499.6 kbp | 2.04 s | 2.85 s | **transition** (gap = 1.4× seed) |
| mm.chr10[40-50M] × rn6.chr20[40-50M] | 3,768.6 kbp | 2.24 s | 15.19 s | gap (7×) |

Same identity tier, same target chromosome, same query chromosome,
same lastz invocation — the bottleneck flips purely as a function of
which 10 Mbp window we sample.

Full step sweep on the densest slice (mm × rn at chr10/chr20 [40-50M]):

| step | Wall | seed_hit | gap_ext | HSPs | bp aligned |
|---:|---:|---:|---:|---:|---:|
| 1   | 47.3 s | 25.6 s **(54 %)** | 21.0 s **(44 %)** | 1767 | 4.15 Mbp |
| 20  | 22.6 s |  3.4 s ( 15 %)    | 18.9 s ( 83 %)    | 1462 | 4.05 Mbp |
| 30  | 21.2 s |  3.0 s ( 14 %)    | 17.8 s ( 84 %)    | 1314 | 4.00 Mbp |
| 50  | 20.0 s |  2.7 s ( 13 %)    | 17.0 s ( 85 %)    | 1180 | 3.92 Mbp |
| 100 | 17.7 s |  2.2 s ( 13 %)    | 15.1 s ( 85 %)    |  991 | 3.77 Mbp |

**Step=1 is itself at the bottleneck transition** (54 % seed / 44 %
gap, essentially balanced), and by step=20 the bottleneck has fully
flipped (15 % / 83 %). This is the cleanest single demonstration we
have that the crossover is a real, narrow boundary.

**Quantitative rule for "which kernel matters"** (per 100 Mbp² input
product, summarized across all five regime points):

| step value | seed = gap density (the threshold) | < threshold | > threshold |
|---|---|---|---|
| step=1   | ~3-4 Mbp aligned   | seed kernel only (sparse synteny, cross-species at scale) | gap kernel only (dense, within-species, primate-primate) |
| step=20  | ~1-2 Mbp aligned   | seed kernel only | gap kernel only |
| step=100 | ~0.3-0.5 Mbp aligned | seed kernel only | gap kernel only |

The threshold is **step-dependent** because seed work shrinks faster
than gap work as `--step` grows. At default `--step=1` the inflection
point sits near 3-4 Mbp; at the more typical cross-species `--step=20`
it sits near 1-2 Mbp; at aggressive `--step=100` it slides down to
~0.3-0.5 Mbp. This means a workload that's seed-dominated at step=1
can become gap-dominated at step=20 simply because we shrunk the
seed budget faster than the gap budget.

**Five-regime overview** (10 Mbp × 10 Mbp slices, sorted by
alignable density):

| Regime | Identity | bp aligned | seed share @ step=20 | gap share @ step=20 |
|---|---:|---:|---:|---:|
| Sparse synteny (hg.chr19 × mm.chr10) | ~80 % blocks | 0.005 Mbp |  87 % | 3 % |
| Cross-species (bird-Z chrZ self) | ~75 % | 0.54 Mbp |  54 % | 42 % |
| **Dense ~85 % (mm.chr10 × rn6.chr20 [40-50M])** | ~85 % | **4.05 Mbp** | **15 %** | **83 %** |
| Within-species (hg.chr1 × hs1.chr1, match15) | ~99.5 % | 11.33 Mbp |  2 % | 94 % |
| Primate-primate (hg.chr19 × rheMac10.chr19) | ~93 % | 31.23 Mbp |  1 % | 99 % |

The list is sorted by density, not identity — and it makes the
threshold visually obvious: bird-Z and "below" are seed-dominated;
mouse-rat dense and "above" are gap-dominated. **Identity does not
order the regimes** (rhesus at 93 % is more gap-dominated than
within-species at 99.5 %, because primate-primate at default cross-
species seed settings overcounts overlapping HSPs).

**Reproduce.** Add `rn6` to the inventory, fetch (~720 MB, ~30 s
download), slice the densest pair, and run the rat sweep:

```bash
# 1) Fetch rn6 (Rnor_6.0).
python3 bench/fetch_genomes.py --only rn6

# 2) Slice the dense-synteny mouse-rat pair (chr10[40,50M] / chr20[40,50M]).
python3 bench/slice_genome.py \
    bench/data_genomes/mm10.chr10.fa \
    /scratch2/shiv1/lastz-bench-data/slices/mm10.chr10_40_50mb.fa \
    --start 40000000 --length 10000000
python3 bench/slice_genome.py \
    bench/data_genomes/rn6.chr20.fa \
    /scratch2/shiv1/lastz-bench-data/slices/rn6.chr20_40_50mb.fa \
    --start 40000000 --length 10000000
ln -sf /scratch2/shiv1/lastz-bench-data/slices/mm10.chr10_40_50mb.fa \
    bench/data_genomes/mm10.chr10_40_50mb.fa
ln -sf /scratch2/shiv1/lastz-bench-data/slices/rn6.chr20_40_50mb.fa \
    bench/data_genomes/rn6.chr20_40_50mb.fa

# 3) Run the sweep (~5 min, single thread).
./step_sweep_rat.sh
```

To reproduce the within-regime density gradient probe (the "same
identity tier, varying density" table above), iterate the same offset
on both species at `--step=100` for a quick smoke result per slice.
Raw artifacts in `bench/results/step-sweep-rat10mb/`.

### Findings worth flagging

1. **Stage-dominance flips with alignment density.** Cross-species at
   moderate divergence: seed_hit_search dominates (91.5% on bird Z
   cross). Self-vs-self: gapped_extension dominates (82.2% on bird Z
   self). Both are "chromosome scale" but they exercise different parts
   of the algorithm. SegAlign's claim is calibrated on the cross
   profile — that's the realistic whole-genome workload.
2. **Crossover happens around 10 Mbp** (for cross-style inputs). At
   small inputs (≤1 Mbp) `gapped_extension` dominates; at 10 Mbp
   seed_hit_search overtakes (72%); at 80 Mbp real bird-cross it owns
   91.5%. This matches SegAlign's headline claim: at chromosome scale,
   the seed-and-filter stage is the right thing to throw on a GPU.
3. **Super-linear scaling of seed_hit_search.** 1 Mbp → 10 Mbp: input
   area grows 100×, total wall grows 25×, but seed_hit_search grows
   66×. Likely cause is the diagonal-hash data structure (re-hashing +
   rehashed-bucket walks scale super-linearly with seed-hit count,
   itself ∝ target_len × query_len for random sequence).
4. **Inner DP loop is 90-97 % of gapped_extension** in every regime
   profiled — confirmed at cycle granularity (rdtsc) across all five
   regimes (sparse mouse-human, bird-Z, mouse-rat, primate, within-
   species). The inner col loop in `ydrop_one_sided_align` runs
   8-17 cyc/cell depending on regime; setup/bounds/trailing/traceback/
   cleanup combined are ≤10 % even in the most extreme case. If we
   GPU-accelerate gapped extension, the inner cell update is the
   only target — everything else is bookkeeping. Within-species is
   the outlier at 8.7 cyc/cell (vs. ~15 elsewhere), likely because
   high identity collapses the inner-loop branches to a tight
   straight-line update. See "Inside `ydrop_one_sided_align` at
   cycle granularity" in §11 for the full per-substage table.
5. **Real DNA is friendlier than random DNA.** 10 Mbp synthetic → 80
   Mbp real bird-cross (60× the seed-search space): seed_hit_search grew
   86×, gapped_extension grew 21×. Soft-masked repeats (~27% of
   galGal6.chrZ) suppress seeding, and surviving HSP density is lower
   in real DNA than in pseudo-random sequence near alignment threshold.
6. **I/O is irrelevant at chr scale.** 7 s out of 73 min (cross),
   14 s out of 239 min (self). No need to switch from FASTA slices to
   direct .2bit reading.
7. **First eu-stack profile (synthetic 1 Mbp × 1 Mbp).** From
   `bench/results/eustack-smoke/`, default lastz, 312 samples in 4.1 s:
   `ydrop_one_sided_align` 57.7%, `find_table_matches` 18.3%,
   `xdrop_extend_seed_hit` 9.3%, `__vdso_gettimeofday` 3.8% (the
   stage-timer overhead is real but small), `add_word` 1.9% (seed table
   build). At 1 Mbp synthetic the gapped extension share is what
   dominates — consistent with the synthetic scaling table above
   (seed_hit_search doesn't take over until ~10 Mbp+). The hottest
   single source line in the seed lookup is `seeds.c:1306` (10.9% of
   all samples) — that's the chain walk inside `find_table_matches`
   (the `for (s = pt->last[word]; s != NULL; s = pt->prev[s-1])` loop).
8. **The bottleneck axis is alignable density, not identity.** Step-
   axis sweeps across five regimes (sparse-synteny mouse-human, bird-Z,
   mouse-rat dense, primate-primate, within-species — see "Step-axis
   Pareto sweeps" and "Pinning the threshold" subsections above) show
   that the seed→gapped crossover is set by how much actually aligns
   per Mbp² of input product, not directly by identity. The mouse-rat
   sweep at ~85 % identity / dense synteny pinned the threshold:
   **at default `--step=20` the seed-share = gap-share crossover sits
   at ~1-2 Mbp aligned per 100 Mbp²**. Below that, the seed kernel
   matters; above that, only `ydrop_one_sided_align` matters. The
   threshold itself is step-dependent — at `--step=1` it shifts up to
   ~3-4 Mbp, at `--step=100` it shifts down to ~0.3-0.5 Mbp. SegAlign's
   "GPU-both-stages" pitch only earns its keep in the narrow density
   band straddling the threshold (bird-Z is the canonical example).

### What we don't have yet

- **A 3-rep variance estimate at chr scale.** The first run hit the
  2-hour SLURM time limit (warmup + measured rep needs ~150 min). The
  sbatch wrapper now defaults to 6h; submit `bird_z_cross` again with
  `WARMUP=0 REPS=3` for a clean variance estimate (~3.5 hours).
- **eu-stack profile of real DNA at 10 Mbp.** ✓ Done — see "Inside
  seed_hit_search" subsection above. The synthetic 1 Mbp baseline is
  also captured at `bench/results/eustack-smoke/`.
- **Cycle-level breakdown of seed_hit_search substages.** ✓ Done via
  `lastz_TS` (the `stage-timing+rdtsc-substages` branch of the lastz
  fork). See "Inside seed_hit_search at cycle granularity" subsection
  above. Raw reports in `bench/results/rdtsc-substages-bird10mb-001/`.
- **Any GPU comparison.** SegAlign is not yet built / run on this
  cluster; that's a follow-up.
- **N-shard CPU baseline.** Standard practice is to shard the target
  by chromosome and run N parallel lastzes. We need that number to be
  fair to LASTZ when later quoting GPU speedups.
- **Dense-array CPU port of the position table.** SegAlign builds a
  flat `(index_table, pos_table)` via counting-sort (see
  `SegAlign/common/seed_pos_table.cu`). Porting just the CPU side of
  that — same algorithm, contiguous instead of linked-list — is the
  obvious first optimization based on the rdtsc finding that the
  chain walk is DRAM-latency bound at 10 Mbp.
- **Gapped-extension cycle breakdown.** ✓ Done — see "Inside
  `ydrop_one_sided_align` at cycle granularity" in §11. The inner DP
  cell loop is 90-97 % of accounted gapped cycles in all five
  regimes at 8-17 cyc/cell; cycle accounting closes to ~70 % of wall,
  consistent with the inner loop being **memory-traffic-bound** more
  than ALU-bound. Raw reports in `bench/results/rdtsc-gapped-smoke/`.
  SegAlign's GPU port handles this kernel via banded x-drop on GPU;
  the next CPU baseline to build is a banded ydrop variant for
  apples-to-apples comparison.
- **rdtsc instrumentation for `find_table_matches_resolve`.** The
  overweight-seed path (any seed with bit-weight > 24, e.g. `match15`)
  bypasses our current substage timers. Small patch to fix; only
  worthwhile if we want substage breakdowns on within-species
  workloads.

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

## 12. Open questions / next steps

In rough order:

1. ✓ **Add `__rdtsc` accumulators to seed_hit_search substages.**
   Done — see the
   [stage-timing+rdtsc-substages branch](https://github.com/shivsundram/lastz/tree/stage-timing+rdtsc-substages)
   and the "Inside seed_hit_search at cycle granularity" subsection of §11.
   The cycle-level data answered both questions raised by the
   eu-stack profile: (i) the chain walk is **memory-latency-bound**
   (89 cyc/iter @ 1 Mbp synth → 125 cyc/iter @ 10 Mbp real, sitting
   between L3 hit and DRAM latency), (ii) the x-drop bar on the
   eu-stack flame is **99.999 % failed extensions**, which is exactly
   the GPU's strong suit (millions of independent identical-shape
   work units, parallelizable across thread blocks).
2. **Prototype a dense-array CPU position table.** SegAlign uses
   counting-sort (`SegAlign/common/seed_pos_table.cu::GenerateSeedPosTable`)
   to convert `pt->last[]` / `pt->prev[]` linked lists into a
   contiguous `(index_table[N_kmers+1], pos_table[N_positions])`
   pair — same algorithm as lastz's chain walk, but reads become
   sequential. Expected to recoup most of the 7 s currently in
   `ftm.chain_walk` on bird-Z 10 Mbp by letting the hardware
   prefetcher hide DRAM latency. Add the `__rdtsc` timers from #1 to
   the new path so we can A/B with apples-to-apples cycle numbers.
3. ✓ **Instrument `gapped_extend.c` with rdtsc substage counters.**
   Done — see "Inside `ydrop_one_sided_align` at cycle granularity"
   in §11. Verdict: the inner DP cell loop is the only GPU kernel
   target; it owns 89.8-96.5 % of accounted gapped cycles across
   all five regimes at 8-17 cyc/cell. Setup/bounds/trailing/
   traceback/cleanup combined are ≤10 % everywhere. The within-
   species regime (hg-vs-CHM13 match15) is 1.9× cheaper per cell
   than the others, suggesting the inner loop simplifies to a
   straight-line update when identity is very high.
3a. ✓ **Two textbook y-drop references (`ydrop_sane`)** both
   validated bit-identical against lastz on 12032 cases (8 hand-
   crafted + 3000 random fuzz pairs, each run forward and reversed,
   each comparison run against both impls). See "Two textbook
   references for `ydrop_one_sided_align`" in §11.
   - `ydrop_one_sided_align_impl_sane`: full `(M+1) × (N+1)`
     matrices, the obvious-but-wasteful version (~52 MB at M=N=2000).
     Reading copy.
   - `ydrop_one_sided_align_impl_sane_double_buffered`: lastz-matching
     memory profile (sweep-row scores + band-compact link tape),
     ~600 KB per call at M=N=2000. Use this for in-situ runtime
     A/B swap against `ydrop_one_sided_align`.
3b. ✓ **In-situ runtime A/B swap** wired into
   `ydrop_one_sided_align` via `YDROP_SANE_IMPL=1`. See "Runtime A/B
   swap inside lastz" in §11. Sane within noise of lastz at 10 Mbp
   (~90 s) on a single contiguous alignment; ~20 % per-call overhead
   amortizes away by then.
3c. ✓ **Cross-anchor sequential-dependency study.** See "Cross-anchor
   sequential-dependency study" in §11. Three orthogonal experiments:
   anchor-loop containment + chain-depth counters, per-call DP-shape
   CSV log (`YDROP_CALL_LOG=path`), and ablation gates
   `LASTZ_DISABLE_NEIGHBOR_MASK=1` / `LASTZ_DISABLE_CONTAINMENT=1`.
   Findings: containment-skip is essential (10× output blowup without
   it, on the same homologous bp set), neighbor masking is cosmetic
   (output bp-set Jaccard 0.9977 vs baseline, +8 % wall when off).
   Median y-drop call is ~500 K cells with a ~600-wide band, far
   below GPU saturation on its own. Path forward: GPU y-drop kernel
   with no masking, persistent threads + work queue to batch calls,
   parallel-greedy-suppress as the GPU post-pass replacement for
   `msp_left_right`.
4. **Re-run bird_z_cross with variance.** Submit again with
   `WARMUP=0 REPS=3` and the new 6h time limit; produces a clean
   median + min/max for the chr-scale wall and stage breakdown.
5. **Run `bird_z_cross_breakdown`** (default + nogapped on the same
   pair) to back out gapped-extension cost on real data, the same way
   we did for synthetic. Should reproduce the 91/8 split.
6. **Build + run SegAlign on the same chrZ pair** for a CPU-vs-GPU
   datapoint. Their docker image avoids the build complexity.
7. **Add an N-shard CPU baseline** to bench.py — launch N parallel
   `lastz` subprocesses on disjoint target slices (chromosome shards),
   produce a "CPU scaling curve" so the GPU comparison isn't 1-core vs
   1-GPU.
8. **Hardware counters.** Either fix `perf` (sysctl
   `kernel.perf_event_paranoid=2`) or rely on the rdtsc-instrumented
   binary above for cycles + an estimate of memory-bandwidth
   utilization. We need to know whether the seed-hit hot loop is ALU-
   or memory-bound before designing a GPU kernel for it.
9. **Whole-genome runs.** Once chr-scale numbers look right, run
   hg38 vs mm10 whole-genome (overnight job). That's the SegAlign
   headline benchmark.

---

## 13. Gotchas, lessons learned, and "don't do that"

A list of things that wasted time during setup so they don't waste
time again:

- **Don't run benchmarks on the head node.** Use `sbatch`. Even
  single-thread numbers will be noisy under shared load. **And `pkill
  -f lastz_T` doesn't reliably kill a Python harness wrapping lastz
  via subprocess — it only kills the lastz child.** The Python parent
  keeps running and respawns. Use `kill <pid>` on the actual python
  PID, or `kill -9 -- -<pgid>` to nuke the whole process group.
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
- **eu-stack can't ptrace siblings under `ptrace_scope=1`.** Even
  though Ubuntu's default policy lets you ptrace your own children, if
  you launch *both* eu-stack and lastz from the same parent, they're
  siblings and eu-stack will fail with `dwfl_thread_getframes:
  Operation not permitted`. The fix is in `profile_eustack.py`: have
  the lastz child call `prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY)` via
  `subprocess.Popen(..., preexec_fn=...)`. Don't refactor that without
  re-testing on a fresh shell.
- **`os.kill(pid, 0)` returns success on zombies.** The first version
  of `profile_eustack.py` used it as a "did the target exit" check and
  ended up sampling a zombie thousands of times after lastz had already
  finished. Use `subprocess.Popen.poll()` instead — that returns the
  exit code as soon as the child becomes a zombie and we can stop the
  sampling loop.
- **`addr2line` on a PIE binary needs file offsets, not vaddrs.**
  `lastz_T` is built as a PIE (default on modern toolchains), so
  eu-stack's reported `0x55....` addresses are virtual addresses with
  the ASLR-randomized base baked in. To resolve to file:line you must
  capture the load base (we read `/proc/<pid>/maps` while the target
  is alive) and subtract it before calling `addr2line`. Symptom of
  forgetting: every line resolves to `??:?`.
- **`eu-stack` doesn't print file:line, only function names.** It's a
  stack walker, not a symbolizer. We do post-resolution via
  `addr2line -e lastz_T -f -C -p <hex offsets>`. If you ever switch
  profilers, anything that already prints file:line per frame will
  Just Work — but eu-stack alone won't.
- **GCC `-O3` partial-inlining mangles file:line attribution.**
  `ydrop_one_sided_align.part.0` (a partial-inlining clone) shows up
  in `top_lines.txt` at `seed_search.c:399` (the call site that pulled
  it in), not `gapped_extend.c` where it's actually defined. Keep this
  in mind when reading the line histogram — function names are still
  correct, only the source location is "the call site".
- **`--extra-arg=<value>` requires the `=` form** for any value that
  starts with `--`. argparse won't accept `--extra-arg --format=maf-`
  as two tokens (it sees `--format=maf-` as another flag). Use
  `--extra-arg=--format=maf-`.

---

## Quick reference card

```bash
# Build
make -C lastz                                                # production
make -C lastz build_lastz_timed                              # lastz_T  (stage timers)
make -C lastz build_lastz_substages                          # lastz_TS (rdtsc substages)

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

# Profile a single run with eu-stack (function + file:line histogram)
python3 bench/profile_eustack.py \
    --lastz lastz/src/lastz_T \
    --target bench/data_genomes/galGal6.chrZ_0_10mb.fa \
    --query  bench/data_genomes/taeGut2.chrZ_0_10mb.fa \
    --out-dir bench/results/eustack-bird10mb-001 \
    --rate-hz 100 --extra-arg=--format=maf-

# Or via SLURM:
sbatch --export=ALL,\
    TARGET=bench/data_genomes/galGal6.chrZ_0_10mb.fa,\
    QUERY=bench/data_genomes/taeGut2.chrZ_0_10mb.fa,\
    OUT_DIR=bench/results/eustack-bird10mb-001 \
    bench/sbatch/profile_eustack.sbatch

# Cycle-level breakdown of seed_hit_search + ydrop_one_sided_align (rdtsc substages)
mkdir -p bench/results/rdtsc-bird10mb-002
LASTZ_STAGE_REPORT=bench/results/rdtsc-bird10mb-002/stage_TS.txt \
    taskset -c 0 lastz/src/lastz_TS \
        bench/data_genomes/galGal6.chrZ_0_10mb.fa \
        bench/data_genomes/taeGut2.chrZ_0_10mb.fa \
        --format=maf- > /dev/null
# Seed substages (chain walk, x-drop, reporter, dedup):
sed -n '/--- rdtsc substage breakdown ---/,/--- rdtsc gapped_extend/p' \
    bench/results/rdtsc-bird10mb-002/stage_TS.txt
# Gapped substages (setup, row bounds, inner DP cells, traceback, cleanup):
sed -n '/--- rdtsc gapped_extend breakdown ---/,/===STAGE_TIMING_END===/p' \
    bench/results/rdtsc-bird10mb-002/stage_TS.txt

# Step-axis Pareto sweeps across five regimes (10 Mbp x 10 Mbp scale)
./step_sweep.sh           # within-species: hg-vs-CHM13, ~75 s
./step_sweep_bird.sh      # cross-species: bird-Z, ~3 min
./step_sweep_mouse.sh     # sparse synteny: mouse-human, ~30 s
./step_sweep_rat.sh       # threshold pin: mouse-rat ~85% dense, ~5 min
./step_sweep_rhesus.sh    # primate-primate (~93%): hg-vs-rhesus, ~17 min
./get_tradeoffs.sh        # within-species headline + bottleneck-flip table

# Inspect
ls   bench/results/
cat  bench/results/<run>/summary.csv | column -ts,
cat  bench/results/<run>/runs/<wl>__<var>__rep0.stage.txt   # raw stage report
cat  bench/results/<run>/top_funcs.txt                      # eu-stack profile
cat  bench/results/<run>/top_lines.txt                      # eu-stack + addr2line
head -10 bench/results/<run>/folded.txt                     # flamegraph input
```

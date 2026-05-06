#!/usr/bin/env python3
"""LASTZ benchmark harness.

Runs (workload x variant) cells under /usr/bin/time -v, collects wall/user/sys
time, peak RSS, page faults, and context switches, and writes structured
results to bench/results/<run-name>/.

Usage examples:

    # List everything
    ./bench.py --list

    # Run the smoke suite, 3 reps each
    ./bench.py --suite smoke --reps 3

    # Run one cell with 5 reps, pinned to CPU 7
    ./bench.py --workload synth.uniform_1000kb --variant default --reps 5 --pin-cpu 7

    # Run a custom name
    ./bench.py --suite scaling --reps 3 --run-name scaling-baseline-2026-05-06

Output layout:

    bench/results/<run_name>/
        manifest.json          # env, lastz version, args, hostname, etc.
        runs/<wl>__<var>__<rep>.json   # per-rep details
        summary.csv            # aggregated table (median across reps)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import re
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DEFAULT_LASTZ = REPO_ROOT / "lastz" / "src" / "lastz"
DEFAULT_LASTZ_TIMED = REPO_ROOT / "lastz" / "src" / "lastz_T"

sys.path.insert(0, str(HERE))
import workloads as wl_mod  # noqa: E402


# ---------------------------------------------------------------------------
# Stage-timing report parsing
# ---------------------------------------------------------------------------

# `lastz_T` (built with -DdbgTiming -DdbgTimingGappedExtend) emits a stage
# timing report bracketed by ===STAGE_TIMING_BEGIN=== / ===STAGE_TIMING_END===
# either to stderr or, if the LASTZ_STAGE_REPORT env var is set, to that file.
# Each line within is "<label>:<ws><value>" where value is either a float
# (seconds) or, for a couple of lines, things like "<int>" or
# "<float> (<float> per second)". We capture all of them and, separately,
# expose the per-stage seconds as a flat dict for easy aggregation.

STAGE_BEGIN_MARKER = "===STAGE_TIMING_BEGIN==="
STAGE_END_MARKER   = "===STAGE_TIMING_END==="

# Map raw label -> snake_case key for the structured "stages" dict. Only the
# numeric per-stage seconds make it here; the per-query and counts go into a
# separate "stage_extras" dict.
_STAGE_LABEL_TO_KEY = {
    "total run time":               "total_run_time_s",
    "sequence 1 I/O":               "sequence_1_io_s",
    "seed position table":          "seed_position_table_s",
    "sequence 2 I/O":               "sequence_2_io_s",
    "seed hit search":              "seed_hit_search_s",
    "chaining":                     "chaining_s",
    "gapped extension":             "gapped_extension_s",
    "interpolation":                "interpolation_s",
    "output":                       "output_s",
    "total query time":             "total_query_time_s",
    # gapped-extension internals (from -DdbgTimingGappedExtend)
    "total time in above_below()":          "ge_above_below_s",
    "total time in left_right()":           "ge_left_right_s",
    "total time in ydrop_align()":          "ge_ydrop_align_s",
    "ydrop_one_sided_align()":              "ge_ydrop_one_sided_align_s",
    "update_lr_bounds()":                   "ge_update_lr_bounds_s",
    "next_sweep_seg()":                     "ge_next_sweep_seg_s",
    "prev_sweep_seg()":                     "ge_prev_sweep_seg_s",
    "update_active_segs()":                 "ge_update_active_segs_s",
    "filter_active_segs()":                 "ge_filter_active_segs_s",
}


def parse_stage_report(text: str) -> tuple[dict, dict]:
    """Parse the stage-timing block, returning (stages_seconds, extras).

    `stages_seconds` maps snake_case keys to floats (seconds). `extras`
    captures everything else verbatim (queries count, per-query rates) so
    nothing is silently dropped.
    """
    stages: dict[str, float] = {}
    extras: dict[str, str] = {}

    if not text:
        return stages, extras

    # Find the bracket; the report may also appear on stderr after marker
    # lines so we just grab the inside of the first occurrence.
    if STAGE_BEGIN_MARKER in text and STAGE_END_MARKER in text:
        start = text.index(STAGE_BEGIN_MARKER) + len(STAGE_BEGIN_MARKER)
        end = text.index(STAGE_END_MARKER, start)
        body = text[start:end]
    else:
        body = text

    for raw in body.splitlines():
        line = raw.rstrip()
        if not line.strip():
            continue
        # Lines look like:  "label:   value"  (label may have parens, etc.)
        # "%-26s %.3f" formatting → split on first ':' with at least one ws.
        m = re.match(r"^\s*(.+?):\s+(.*)$", line)
        if not m:
            continue
        label, value = m.group(1).strip(), m.group(2).strip()
        if label in _STAGE_LABEL_TO_KEY:
            try:
                stages[_STAGE_LABEL_TO_KEY[label]] = float(value)
            except ValueError:
                extras[label] = value
        else:
            extras[label] = value

    return stages, extras


def read_stage_report_file(path: Path) -> tuple[dict, dict]:
    if not path.exists():
        return {}, {}
    try:
        return parse_stage_report(path.read_text())
    except Exception as e:  # malformed report shouldn't kill the whole run
        return {}, {"_parse_error": str(e)}


# ---------------------------------------------------------------------------
# /usr/bin/time -v parsing
# ---------------------------------------------------------------------------

# `/usr/bin/time -v` writes a fixed set of "Key: value" lines to stderr.
# Example excerpt:
#   Command being timed: "lastz a b"
#   User time (seconds): 0.12
#   System time (seconds): 0.00
#   Percent of CPU this job got: 99%
#   Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.13
#   Maximum resident set size (kbytes): 2872
#   Voluntary context switches: 1
#   Involuntary context switches: 3
#   Major (requiring I/O) page faults: 0
#   Minor (reclaiming a frame) page faults: 145

_TIME_KEY_MAP = {
    "User time (seconds)":              ("user_s",            float),
    "System time (seconds)":            ("sys_s",             float),
    "Percent of CPU this job got":      ("cpu_pct",           lambda s: int(s.rstrip("%"))),
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": ("elapsed_str", str),
    "Maximum resident set size (kbytes)": ("max_rss_kb",      int),
    "Voluntary context switches":       ("ctxsw_voluntary",   int),
    "Involuntary context switches":     ("ctxsw_involuntary", int),
    "Major (requiring I/O) page faults": ("page_faults_major", int),
    "Minor (reclaiming a frame) page faults": ("page_faults_minor", int),
    "Exit status":                      ("time_reported_exit", int),
}


def parse_elapsed(elapsed_str: str) -> float:
    """Parse `/usr/bin/time -v`'s elapsed field into seconds."""
    parts = elapsed_str.strip().split(":")
    if len(parts) == 3:
        h, m, s = parts
        return int(h) * 3600 + int(m) * 60 + float(s)
    if len(parts) == 2:
        m, s = parts
        return int(m) * 60 + float(s)
    return float(parts[0])


def parse_time_v(stderr_text: str) -> dict:
    # Some keys in `/usr/bin/time -v` contain ':' characters (e.g.
    # "Elapsed (wall clock) time (h:mm:ss or m:ss)"), so partitioning on the
    # first ':' is wrong. Match by longest-known-key prefix instead.
    out: dict = {}
    sorted_keys = sorted(_TIME_KEY_MAP.keys(), key=len, reverse=True)
    for raw in stderr_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        for key in sorted_keys:
            prefix = key + ":"
            if line.startswith(prefix):
                val = line[len(prefix):].strip()
                field_name, caster = _TIME_KEY_MAP[key]
                try:
                    out[field_name] = caster(val)
                except (ValueError, TypeError):
                    out[field_name] = val
                break
    if "elapsed_str" in out:
        out["wall_s"] = parse_elapsed(out["elapsed_str"])
    return out


# ---------------------------------------------------------------------------
# Run model
# ---------------------------------------------------------------------------

@dataclass
class RunResult:
    workload: str
    variant: str
    rep: int
    cmd: list[str]
    target: str
    query: str
    extra_args: list[str]
    exit_code: int
    timed_out: bool
    wall_s: float | None = None
    user_s: float | None = None
    sys_s: float | None = None
    cpu_pct: int | None = None
    max_rss_kb: int | None = None
    page_faults_major: int | None = None
    page_faults_minor: int | None = None
    ctxsw_voluntary: int | None = None
    ctxsw_involuntary: int | None = None
    output_bytes: int = 0
    stderr_tail: str = ""
    timestamp: str = ""
    workload_meta: dict = field(default_factory=dict)
    # Stage-timing report (only populated when running `lastz_T`). `stages`
    # holds floats in seconds; `stage_extras` keeps the per-query counts and
    # any unrecognized lines verbatim.
    stages: dict = field(default_factory=dict)
    stage_extras: dict = field(default_factory=dict)


def discover_lastz(explicit: Path | None, prefer_timed: bool) -> tuple[Path, bool]:
    """Locate a lastz binary. Returns (path, is_instrumented)."""
    if explicit is not None:
        if not explicit.exists():
            raise FileNotFoundError(f"lastz binary not found: {explicit}")
        # Heuristic: any binary at .../lastz_T is treated as instrumented.
        return explicit, explicit.name.endswith("_T")
    if prefer_timed and DEFAULT_LASTZ_TIMED.exists():
        return DEFAULT_LASTZ_TIMED, True
    if DEFAULT_LASTZ.exists():
        return DEFAULT_LASTZ, False
    if DEFAULT_LASTZ_TIMED.exists():
        return DEFAULT_LASTZ_TIMED, True
    on_path = shutil.which("lastz")
    if on_path:
        return Path(on_path), False
    raise FileNotFoundError(
        f"Could not find lastz. Tried {DEFAULT_LASTZ}, {DEFAULT_LASTZ_TIMED}, and $PATH."
        " Build it (cd lastz && make build_lastz_timed) or pass --lastz."
    )


def lastz_version_string(lastz: Path) -> str:
    try:
        r = subprocess.run([str(lastz), "--version"],
                           capture_output=True, text=True, timeout=10)
        return (r.stdout or r.stderr).strip()
    except Exception as e:
        return f"<unavailable: {e}>"


def _spec(path: Path, actions: Iterable[str]) -> str:
    """Build a lastz input spec, e.g. '/path/to/x.fa[multiple,subset=foo]'."""
    actions = list(actions)
    if not actions:
        return str(path)
    return f"{path}[{','.join(actions)}]"


def build_command(
    lastz: Path,
    target_spec: str,
    query_spec: str,
    workload_args: Iterable[str],
    variant_args: Iterable[str],
    output_path: Path,
    pin_cpu: int | None,
) -> list[str]:
    cmd: list[str] = []
    if pin_cpu is not None:
        cmd += ["taskset", "-c", str(pin_cpu)]
    cmd += ["/usr/bin/time", "-v",
            str(lastz), target_spec, query_spec,
            f"--output={output_path}"]
    cmd += list(workload_args)
    cmd += list(variant_args)
    return cmd


def run_one(
    lastz: Path,
    instrumented: bool,
    workload: wl_mod.Workload,
    variant: wl_mod.Variant,
    rep: int,
    out_dir: Path,
    pin_cpu: int | None,
    timeout_s: float | None,
    keep_output: bool,
) -> RunResult:
    target, query, wl_meta = workload.resolve()
    target_spec = _spec(target, workload.target_actions)
    query_spec  = _spec(query,  workload.query_actions)

    # Per-rep output file under runs/ alongside the JSON.
    output_path = out_dir / f"{workload.id}__{variant.id}__{rep}.lastz_out"
    cmd = build_command(lastz, target_spec, query_spec,
                        workload.lastz_args, variant.extra_args,
                        output_path, pin_cpu)

    timestamp = time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime())

    # If instrumented, route the stage report into a per-rep file via env var
    # so it doesn't tangle with regular stderr (and so we don't have to scrape
    # multi-line markers out of stderr).
    env = os.environ.copy()
    stage_path: Path | None = None
    if instrumented:
        stage_path = out_dir / f"{workload.id}__{variant.id}__{rep}.stage.txt"
        env["LASTZ_STAGE_REPORT"] = str(stage_path)

    timed_out = False
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            env=env,
        )
        exit_code = proc.returncode
        stderr_text = proc.stderr
    except subprocess.TimeoutExpired as e:
        timed_out = True
        exit_code = -1
        stderr_text = (e.stderr or "") if isinstance(e.stderr, str) else ""

    parsed = parse_time_v(stderr_text)

    stages: dict = {}
    stage_extras: dict = {}
    if instrumented and stage_path is not None:
        stages, stage_extras = read_stage_report_file(stage_path)

    output_bytes = output_path.stat().st_size if output_path.exists() else 0

    # Trim stderr — drop the time -v block to leave just lastz's own diagnostics.
    lastz_stderr_lines = [
        ln for ln in stderr_text.splitlines()
        if not any(ln.strip().startswith(k) for k in _TIME_KEY_MAP)
        and not ln.strip().startswith("Command being timed:")
        and not ln.strip().startswith("Average")
        and not ln.strip().startswith("File system")
        and not ln.strip().startswith("Socket")
        and not ln.strip().startswith("Signals")
        and not ln.strip().startswith("Page size")
        and not ln.strip().startswith("Swaps")
        and not ln.strip().startswith("Maximum")  # caught above too, just in case
    ]
    stderr_tail = "\n".join(lastz_stderr_lines[-20:])

    if not keep_output and output_path.exists():
        output_path.unlink()

    return RunResult(
        workload=workload.id,
        variant=variant.id,
        rep=rep,
        cmd=cmd,
        target=str(target),
        query=str(query),
        extra_args=list(variant.extra_args),
        exit_code=exit_code,
        timed_out=timed_out,
        wall_s=parsed.get("wall_s"),
        user_s=parsed.get("user_s"),
        sys_s=parsed.get("sys_s"),
        cpu_pct=parsed.get("cpu_pct"),
        max_rss_kb=parsed.get("max_rss_kb"),
        page_faults_major=parsed.get("page_faults_major"),
        page_faults_minor=parsed.get("page_faults_minor"),
        ctxsw_voluntary=parsed.get("ctxsw_voluntary"),
        ctxsw_involuntary=parsed.get("ctxsw_involuntary"),
        output_bytes=output_bytes,
        stderr_tail=stderr_tail,
        timestamp=timestamp,
        workload_meta=wl_meta,
        stages=stages,
        stage_extras=stage_extras,
    )


# ---------------------------------------------------------------------------
# Aggregation + CSV
# ---------------------------------------------------------------------------

BASE_SUMMARY_FIELDS = [
    "workload", "variant", "n_reps", "ok",
    "wall_s_min", "wall_s_median", "wall_s_max",
    "user_s_median", "sys_s_median",
    "cpu_pct_median",
    "max_rss_kb_median",
    "page_faults_major_median", "page_faults_minor_median",
    "output_bytes_median",
]

# Per-stage columns are appended dynamically based on which stage keys appear
# in the run results (so a non-instrumented run still produces a usable CSV).
STAGE_KEY_ORDER = [
    "total_run_time_s",
    "sequence_1_io_s",
    "seed_position_table_s",
    "sequence_2_io_s",
    "seed_hit_search_s",
    "chaining_s",
    "gapped_extension_s",
    "interpolation_s",
    "output_s",
    "total_query_time_s",
    "ge_above_below_s",
    "ge_left_right_s",
    "ge_ydrop_align_s",
    "ge_ydrop_one_sided_align_s",
    "ge_update_lr_bounds_s",
    "ge_next_sweep_seg_s",
    "ge_prev_sweep_seg_s",
    "ge_update_active_segs_s",
    "ge_filter_active_segs_s",
]


def _median_or_none(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return None
    return statistics.median(xs)


def summarise(results: list[RunResult]) -> tuple[list[dict], list[str]]:
    """Aggregate per-cell medians. Returns (rows, fieldnames)."""
    by_cell: dict[tuple[str, str], list[RunResult]] = {}
    for r in results:
        by_cell.setdefault((r.workload, r.variant), []).append(r)

    # Determine which stage keys appear at least once so we only emit columns
    # that have any data. Preserves canonical order.
    seen_stage_keys = set()
    for r in results:
        seen_stage_keys.update(r.stages.keys())
    stage_cols = [k for k in STAGE_KEY_ORDER if k in seen_stage_keys]

    rows: list[dict] = []
    for (wl, var), rs in by_cell.items():
        wall = [r.wall_s for r in rs if r.wall_s is not None]
        ok = all(r.exit_code == 0 and not r.timed_out for r in rs)
        row = {
            "workload": wl,
            "variant": var,
            "n_reps": len(rs),
            "ok": int(ok),
            "wall_s_min": min(wall) if wall else None,
            "wall_s_median": statistics.median(wall) if wall else None,
            "wall_s_max": max(wall) if wall else None,
            "user_s_median": _median_or_none([r.user_s for r in rs]),
            "sys_s_median":  _median_or_none([r.sys_s for r in rs]),
            "cpu_pct_median": _median_or_none([r.cpu_pct for r in rs]),
            "max_rss_kb_median": _median_or_none([r.max_rss_kb for r in rs]),
            "page_faults_major_median": _median_or_none([r.page_faults_major for r in rs]),
            "page_faults_minor_median": _median_or_none([r.page_faults_minor for r in rs]),
            "output_bytes_median": _median_or_none([r.output_bytes for r in rs]),
        }
        for k in stage_cols:
            row[f"{k}_median"] = _median_or_none([r.stages.get(k) for r in rs])
        rows.append(row)
    rows.sort(key=lambda r: (r["workload"], r["variant"]))
    fieldnames = BASE_SUMMARY_FIELDS + [f"{k}_median" for k in stage_cols]
    return rows, fieldnames


def write_summary_csv(rows: list[dict], fieldnames: list[str], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_list() -> int:
    print("Workloads:")
    for wid, w in sorted(wl_mod.WORKLOADS.items()):
        tags = f" [{', '.join(w.tags)}]" if w.tags else ""
        print(f"  {wid}{tags}")
        print(f"      {w.description}")
    print("\nVariants:")
    for vid, v in sorted(wl_mod.VARIANTS.items()):
        args = " ".join(v.extra_args) if v.extra_args else "(no extra args)"
        print(f"  {vid}: {v.description}")
        print(f"      args: {args}")
    print("\nSuites:")
    for sid, cells in sorted(wl_mod.SUITES.items()):
        print(f"  {sid}  ({len(cells)} cell{'s' if len(cells) != 1 else ''})")
        for wlid, varid in cells:
            print(f"      {wlid}  +  {varid}")
    return 0


def collect_cells(args: argparse.Namespace) -> list[tuple[str, str]]:
    cells: list[tuple[str, str]] = []
    if args.suite:
        if args.suite not in wl_mod.SUITES:
            raise SystemExit(f"unknown suite: {args.suite}")
        cells.extend(wl_mod.SUITES[args.suite])
    if args.workload:
        var = args.variant or "default"
        cells.append((args.workload, var))
    if not cells:
        raise SystemExit("nothing to run; pass --suite or --workload")
    # Validate
    for wlid, varid in cells:
        if wlid not in wl_mod.WORKLOADS:
            raise SystemExit(f"unknown workload: {wlid}")
        if varid not in wl_mod.VARIANTS:
            raise SystemExit(f"unknown variant: {varid}")
    return cells


def env_manifest(lastz: Path, instrumented: bool, args: argparse.Namespace) -> dict:
    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cpu_count": os.cpu_count(),
        "lastz_path": str(lastz),
        "lastz_instrumented": instrumented,
        "lastz_version": lastz_version_string(lastz),
        "argv": sys.argv,
        "reps": args.reps,
        "warmup": args.warmup,
        "pin_cpu": args.pin_cpu,
        "timeout_s": args.timeout_s,
        "keep_output": args.keep_output,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime()),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list", action="store_true", help="List workloads/variants/suites and exit.")
    p.add_argument("--suite", type=str, default=None)
    p.add_argument("--workload", type=str, default=None)
    p.add_argument("--variant", type=str, default=None,
                   help="Defaults to 'default' when --workload is given.")
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--warmup", type=int, default=1, help="Warm-up reps (results discarded).")
    p.add_argument("--pin-cpu", type=int, default=None,
                   help="If set, pin lastz to this CPU id via taskset.")
    p.add_argument("--lastz", type=Path, default=None,
                   help=f"Path to lastz binary (default: prefer {DEFAULT_LASTZ_TIMED}, "
                        f"fall back to {DEFAULT_LASTZ})")
    p.add_argument("--no-stage-timers", action="store_true",
                   help="Use the un-instrumented lastz binary even if lastz_T is available.")
    p.add_argument("--results-dir", type=Path, default=HERE / "results")
    p.add_argument("--run-name", type=str, default=None,
                   help="Subdirectory under results/. Defaults to a timestamp.")
    p.add_argument("--timeout-s", type=float, default=None,
                   help="Per-run wall-clock timeout (seconds).")
    p.add_argument("--keep-output", action="store_true",
                   help="Don't delete lastz output files after each run.")
    args = p.parse_args(argv)

    if args.list:
        return cmd_list()

    lastz, instrumented = discover_lastz(args.lastz, prefer_timed=not args.no_stage_timers)
    cells = collect_cells(args)

    run_name = args.run_name or time.strftime("run-%Y%m%d-%H%M%S")
    run_dir = args.results_dir / run_name
    runs_dir = run_dir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    manifest = env_manifest(lastz, instrumented, args)
    manifest["cells"] = [{"workload": w, "variant": v} for w, v in cells]
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"# lastz: {lastz}  (stage_timers={'on' if instrumented else 'off'})")
    print(f"# version: {manifest['lastz_version'].splitlines()[0]}")
    print(f"# results: {run_dir}")
    print(f"# cells: {len(cells)}, reps: {args.reps} (+{args.warmup} warmup)")
    print()

    all_results: list[RunResult] = []
    for wlid, varid in cells:
        workload = wl_mod.WORKLOADS[wlid]
        variant = wl_mod.VARIANTS[varid]

        # Warm-up reps (discarded) — useful for filling page cache so the
        # measured reps see steady-state I/O behavior.
        for w in range(args.warmup):
            _ = run_one(lastz, instrumented, workload, variant, rep=-(w + 1),
                        out_dir=runs_dir, pin_cpu=args.pin_cpu,
                        timeout_s=args.timeout_s, keep_output=False)

        cell_results: list[RunResult] = []
        for rep in range(args.reps):
            r = run_one(lastz, instrumented, workload, variant, rep=rep,
                        out_dir=runs_dir, pin_cpu=args.pin_cpu,
                        timeout_s=args.timeout_s,
                        keep_output=args.keep_output)
            cell_results.append(r)
            json_path = runs_dir / f"{wlid}__{varid}__{rep}.json"
            json_path.write_text(json.dumps(asdict(r), indent=2, default=str))

            status = "OK" if r.exit_code == 0 and not r.timed_out else f"FAIL exit={r.exit_code}"
            ge = r.stages.get("gapped_extension_s")
            ge_frac = (ge / r.wall_s) if (ge is not None and r.wall_s) else None
            ge_str = f"  ge={ge:>6.3f}s ({ge_frac*100:4.1f}%)" if ge_frac is not None else ""
            print(f"  {wlid:40s}  {varid:14s}  rep={rep}  "
                  f"wall={r.wall_s if r.wall_s is not None else 'NA':>7}s  "
                  f"rss={r.max_rss_kb if r.max_rss_kb is not None else 'NA':>8}kb  "
                  f"out={r.output_bytes:>9}B{ge_str}  {status}")
            if r.exit_code != 0 and r.stderr_tail:
                print("    stderr tail:")
                for ln in r.stderr_tail.splitlines()[-5:]:
                    print(f"      {ln}")
        all_results.extend(cell_results)

    summary, fieldnames = summarise(all_results)
    write_summary_csv(summary, fieldnames, run_dir / "summary.csv")

    print()
    print("# summary (medians across reps):")
    if instrumented:
        hdr = (f"  {'workload':40s}  {'variant':12s}  {'n':>3s}  "
               f"{'wall':>7s}  {'seed_idx':>8s}  {'seed_hit':>8s}  "
               f"{'gapped':>8s}  {'output':>7s}  {'ydrop1s':>8s}")
    else:
        hdr = (f"  {'workload':40s}  {'variant':12s}  {'n':>3s}  "
               f"{'wall':>7s}  {'rss_kb':>8s}  {'out_B':>9s}")
    print(hdr)
    for row in summary:
        wm = row["wall_s_median"]
        n = row["n_reps"]
        if instrumented:
            def _f(k):  # `k` should include the `_s` suffix
                v = row.get(f"{k}_median")
                return f"{v:.3f}" if v is not None else "  NA"
            print(f"  {row['workload']:40s}  {row['variant']:12s}  {n:>3d}  "
                  f"{(f'{wm:.3f}' if wm is not None else '  NA'):>7s}  "
                  f"{_f('seed_position_table_s'):>8s}  {_f('seed_hit_search_s'):>8s}  "
                  f"{_f('gapped_extension_s'):>8s}  {_f('output_s'):>7s}  "
                  f"{_f('ge_ydrop_one_sided_align_s'):>8s}")
        else:
            rm = row["max_rss_kb_median"]
            om = row["output_bytes_median"]
            print(f"  {row['workload']:40s}  {row['variant']:12s}  {n:>3d}  "
                  f"{(f'{wm:.3f}' if wm is not None else '  NA'):>7s}  "
                  f"{(f'{int(rm)}' if rm is not None else '  NA'):>8s}  "
                  f"{(f'{int(om)}' if om is not None else '  NA'):>9s}")
    print()
    print(f"# wrote {run_dir / 'summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

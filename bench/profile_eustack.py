#!/usr/bin/env python3
"""eu-stack-based PC sampling profiler for lastz.

Launches lastz_T as a subprocess (so we own its ptrace permission under the
default Ubuntu ptrace_scope=1) and samples its call stack with eu-stack at
a fixed rate. After the target exits, aggregates samples into a folded
stack-trace format compatible with Brendan Gregg's flamegraph.pl.

Output directory layout:

    out_dir/
        manifest.json              # what we ran, when, how many samples
        lastz.out                  # lastz alignment output (kept for sanity)
        lastz.err                  # lastz stderr (incl. stage-timing report)
        raw/sample-NNNNN.txt       # each eu-stack snapshot
        folded.txt                 # frame1;frame2;... <count>  (one line per
                                   # unique stack)
        top_funcs.txt              # top functions by sample count
        top_lines.txt              # top (func, file:line) by sample count

Usage:

    python3 bench/profile_eustack.py \\
        --lastz lastz/src/lastz_T \\
        --target bench/data_genomes/galGal6.chrZ_0_10mb.fa \\
        --query  bench/data_genomes/taeGut2.chrZ_0_10mb.fa \\
        --out-dir bench/results/eustack-bird10mb-001/ \\
        --rate-hz 100 \\
        --extra-arg --format=maf-

Notes:
- This is a sampler, not a tracer; the target pays only the cost of being
  ptrace-attached and walked once per sample (~few ms each), independent of
  workload size.
- We tolerate eu-stack failing on a given sample (process exited mid-call,
  transient ptrace race, etc.) — those samples are dropped.
- Symbol resolution requires lastz_T to be built with -g (it is by default
  in our Makefile patch).
"""

from __future__ import annotations

import argparse
import collections
import ctypes
import ctypes.util
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

# prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY, ...) — declare that ANY process
# may ptrace us. Needed because under Ubuntu's default ptrace_scope=1 only
# direct ancestors can attach; eu-stack would be a sibling of lastz (both
# children of this script), so we have to opt in explicitly.
_PR_SET_PTRACER = 0x59616D61  # "Yama" magic from linux/prctl.h
_PR_SET_PTRACER_ANY = ctypes.c_ulong(-1).value  # 0xFFFFFFFFFFFFFFFF
_LIBC = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)


def _allow_any_ptracer() -> None:
    """preexec_fn for the lastz child: opt in to being ptrace'd by anyone."""
    _LIBC.prctl(_PR_SET_PTRACER, ctypes.c_ulong(_PR_SET_PTRACER_ANY), 0, 0, 0)


# eu-stack output on this system looks like:
#   PID 12345 - process
#   TID 12345:
#   #0  0x00005639edb2e1a2 find_table_matches
#   #1  0x00005639edb2f02c seed_hit_search
#   ...
# (No file:line info — eu-stack is a stack walker, not addr2line. We capture
# the hex PC and resolve to file:line in a separate post-processing pass.)
FRAME_RE = re.compile(
    r"^\s*#\d+\s+0x(?P<addr>[0-9a-fA-F]+)\s+(?P<func>\S+)"
)


def parse_eustack_sample(text: str) -> list[tuple[str, str]] | None:
    """Returns innermost-first list of (func, hex_addr) frames, or None if
    the sample contains no usable stack."""
    frames: list[tuple[str, str]] = []
    for line in text.splitlines():
        m = FRAME_RE.match(line)
        if not m:
            continue
        frames.append((m.group("func"), m.group("addr")))
    return frames if frames else None


def fold_stack(frames: list[tuple[str, str]]) -> str:
    """Innermost-first input → outermost-first semicolon-joined string."""
    return ";".join(f for f, _ in reversed(frames))


def capture_load_base(pid: int, binary: Path) -> int | None:
    """Read /proc/<pid>/maps and find the load base of `binary`.

    For a PIE binary the executable is mapped at an ASLR-randomized address;
    eu-stack reports virtual addresses, but addr2line needs file offsets
    (vaddr - load_base). We grab that base while the target is still alive.
    """
    try:
        text = Path(f"/proc/{pid}/maps").read_text()
    except OSError:
        return None
    # Find the FIRST entry that maps `binary` (highest precedence: text segment).
    binary_str = str(binary.resolve())
    for line in text.splitlines():
        parts = line.split(maxsplit=5)
        if len(parts) >= 6 and parts[5] == binary_str and "x" in parts[1]:
            start = parts[0].split("-")[0]
            return int(start, 16)
    return None


def resolve_addrs(addrs: list[str], binary: Path,
                  load_base: int | None) -> dict[str, str]:
    """Run addr2line once with all unique addresses to get file:line per addr.
    Best-effort: returns empty mapping if addr2line fails or the binary lacks
    DWARF info. For PIE binaries, subtract load_base to get file offsets."""
    if not addrs:
        return {}
    if subprocess.run(["which", "addr2line"], capture_output=True).returncode != 0:
        return {}
    # Convert vaddrs to file offsets if we know the PIE load base.
    if load_base is not None:
        adjusted = [f"0x{int(a, 16) - load_base:x}" for a in addrs]
    else:
        adjusted = [f"0x{a}" for a in addrs]
    try:
        cp = subprocess.run(
            ["addr2line", "-e", str(binary), "-f", "-C", "-p"] + adjusted,
            capture_output=True, text=True, timeout=60.0,
        )
    except subprocess.TimeoutExpired:
        return {}
    if cp.returncode != 0:
        return {}
    out: dict[str, str] = {}
    for addr, line in zip(addrs, cp.stdout.splitlines()):
        # Format with -f -p: "func at file.c:123" OR "func at ??:?"
        m = re.search(r"\sat\s+(?P<loc>.+)$", line)
        loc = m.group("loc").strip() if m else "?"
        if loc.startswith("?"):
            loc = ""
        out[addr] = loc
    return out


def aggregate(raw_dir: Path, out_dir: Path, binary: Path | None,
              load_base: int | None) -> dict:
    folded_counts: collections.Counter[str] = collections.Counter()
    func_counts: collections.Counter[str] = collections.Counter()
    leaf_counts: collections.Counter[tuple[str, str]] = collections.Counter()
    all_leaf_addrs: set[str] = set()
    n_samples = 0
    n_failed = 0

    for path in sorted(raw_dir.glob("sample-*.txt")):
        text = path.read_text(errors="replace")
        frames = parse_eustack_sample(text)
        if not frames:
            n_failed += 1
            continue
        n_samples += 1
        folded_counts[fold_stack(frames)] += 1
        leaf_func, leaf_addr = frames[0]
        func_counts[leaf_func] += 1
        leaf_counts[(leaf_func, leaf_addr)] += 1
        all_leaf_addrs.add(leaf_addr)

    addr_to_loc = (resolve_addrs(sorted(all_leaf_addrs), binary, load_base)
                   if binary else {})

    (out_dir / "folded.txt").write_text(
        "\n".join(f"{stack} {count}" for stack, count in folded_counts.most_common())
        + "\n"
    )
    with open(out_dir / "top_funcs.txt", "w") as f:
        f.write(f"# total samples: {n_samples}  (eu-stack failures: {n_failed})\n")
        f.write(f"{'count':>8} {'pct':>6}  function\n")
        for func, c in func_counts.most_common(50):
            pct = 100.0 * c / n_samples if n_samples else 0
            f.write(f"{c:>8} {pct:>5.1f}%  {func}\n")
    with open(out_dir / "top_lines.txt", "w") as f:
        any_loc = any(addr_to_loc.get(a, "") for _, a in leaf_counts)
        f.write(f"# total samples: {n_samples}  (eu-stack failures: {n_failed})\n")
        if not any_loc:
            f.write("# NOTE: no line info resolved (addr2line failed or binary "
                    "lacks DWARF). Falling back to (func, hex address).\n")
        f.write(f"{'count':>8} {'pct':>6}  function  (location/addr)\n")
        for (func, addr), c in leaf_counts.most_common(50):
            pct = 100.0 * c / n_samples if n_samples else 0
            loc = addr_to_loc.get(addr, "") or f"0x{addr}"
            f.write(f"{c:>8} {pct:>5.1f}%  {func}  ({loc})\n")
    return {"n_samples": n_samples, "n_failed": n_failed,
            "n_unique_stacks": len(folded_counts)}


def sample_loop(proc: subprocess.Popen, out_raw: Path, interval_s: float,
                max_samples: int | None) -> int:
    """Loop sampling proc.pid with eu-stack until the target exits or we hit
    max_samples. Returns the count of attempted samples.

    Uses proc.poll() (not os.kill) so we correctly detect zombies — once
    lastz exits, poll() returns its exit code and we stop sampling.
    """
    pid = proc.pid
    n = 0
    while True:
        if max_samples is not None and n >= max_samples:
            return n
        if proc.poll() is not None:
            return n
        sample_path = out_raw / f"sample-{n:06d}.txt"
        with open(sample_path, "wb") as f:
            try:
                subprocess.run(["eu-stack", "-p", str(pid)],
                               stdout=f, stderr=subprocess.DEVNULL,
                               timeout=2.0)
            except subprocess.TimeoutExpired:
                f.close()
                sample_path.unlink(missing_ok=True)
        n += 1
        time.sleep(interval_s)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lastz", type=Path, required=True,
                   help="Path to the lastz binary (likely lastz_T).")
    p.add_argument("--target", type=Path, required=True)
    p.add_argument("--query",  type=Path, required=True)
    p.add_argument("--extra-arg", action="append", default=[],
                   help="Extra arg passed to lastz (repeat for each).")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--rate-hz", type=float, default=100.0,
                   help="Sample rate in Hz (default: 100).")
    p.add_argument("--max-samples", type=int, default=None,
                   help="Stop sampling after this many attempts (default: until "
                        "the target exits).")
    p.add_argument("--no-keep-raw", action="store_true",
                   help="Delete the per-sample raw eu-stack files after "
                        "aggregation (default: keep — they're small and "
                        "useful for re-aggregating with different filters).")
    args = p.parse_args(argv)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    # Sanity: eu-stack must be available
    if subprocess.run(["which", "eu-stack"], capture_output=True).returncode != 0:
        print("ERROR: eu-stack not found in PATH", file=sys.stderr)
        return 2

    # We launch lastz_T as our subprocess so ptrace_scope=1 lets us sample it.
    cmd = [str(args.lastz), str(args.target), str(args.query)] + args.extra_arg
    started_at = time.time()
    print(f"# launching: {' '.join(cmd)}", file=sys.stderr)
    print(f"# out_dir:   {out_dir}", file=sys.stderr)
    print(f"# rate:      {args.rate_hz} Hz "
          f"(interval {1000.0 / args.rate_hz:.1f} ms)", file=sys.stderr)

    out_log = open(out_dir / "lastz.out", "wb")
    err_log = open(out_dir / "lastz.err", "wb")
    # preexec_fn runs in the child between fork and exec — perfect spot for
    # the PR_SET_PTRACER opt-in.
    proc = subprocess.Popen(cmd, stdout=out_log, stderr=err_log,
                            preexec_fn=_allow_any_ptracer)
    interval_s = 1.0 / args.rate_hz

    # Race-tolerant: a brief settling delay lets lastz set up before we sample.
    time.sleep(0.05)

    # Capture the PIE load base while the process is alive — we need it later
    # to convert eu-stack's vaddrs to file offsets for addr2line.
    load_base = capture_load_base(proc.pid, args.lastz)
    if load_base is not None:
        print(f"# load base: 0x{load_base:x} (lastz_T text segment)", file=sys.stderr)
    else:
        print("# warning: could not capture load base; line resolution will fail",
              file=sys.stderr)

    try:
        n_attempted = sample_loop(proc, raw_dir, interval_s, args.max_samples)
    except KeyboardInterrupt:
        proc.terminate()
        n_attempted = -1

    rc = proc.wait()
    out_log.close()
    err_log.close()
    duration = time.time() - started_at

    print(f"# lastz exited rc={rc} after {duration:.1f}s", file=sys.stderr)
    print(f"# attempted {n_attempted} samples", file=sys.stderr)

    stats = aggregate(raw_dir, out_dir, args.lastz, load_base)

    manifest = {
        "lastz": str(args.lastz),
        "target": str(args.target),
        "query":  str(args.query),
        "extra_args": args.extra_arg,
        "rate_hz": args.rate_hz,
        "interval_ms": 1000.0 / args.rate_hz,
        "max_samples": args.max_samples,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(started_at)),
        "wall_s": duration,
        "lastz_exit": rc,
        "samples_attempted": n_attempted,
        "samples_with_stack": stats["n_samples"],
        "samples_failed": stats["n_failed"],
        "unique_stacks": stats["n_unique_stacks"],
        "load_base_hex": f"0x{load_base:x}" if load_base is not None else None,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # Echo the headline table to stdout for human-readable feedback.
    print()
    print("==== eu-stack profile ====")
    print(f"target  : {args.target}")
    print(f"query   : {args.query}")
    print(f"wall    : {duration:.1f} s")
    print(f"samples : {stats['n_samples']} usable "
          f"({stats['n_failed']} failed) at {args.rate_hz} Hz")
    print(f"unique  : {stats['n_unique_stacks']} distinct call stacks")
    print()
    print("--- top 15 functions (innermost frame at sample time) ---")
    text = (out_dir / "top_funcs.txt").read_text().splitlines()[:17]
    for line in text:
        print(line)

    if args.no_keep_raw:
        for p in raw_dir.glob("sample-*.txt"):
            p.unlink()
        try:
            raw_dir.rmdir()
        except OSError:
            pass

    return 0 if rc == 0 else rc


if __name__ == "__main__":
    raise SystemExit(main())

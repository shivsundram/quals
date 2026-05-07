"""Workload + variant registry for the lastz benchmark harness.

A *workload* defines an input pair (target + query). A *variant* defines a
set of lastz CLI flags applied on top. We measure (workload x variant) cells.

Workloads come in two flavors:

  * "static": a pair of FASTA/2bit files that already exist on disk
    (e.g. the ones in lastz/test_data).
  * "synthetic": generated on demand by gen_data.make_pair, parameterised
    by target/query length and divergence. These let us scale the input
    and watch how runtime + memory behave as a function of size.

Suites are just named lists of workload IDs, for convenience.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

REPO_ROOT = Path(__file__).resolve().parent.parent
LASTZ_TEST_DATA = REPO_ROOT / "lastz" / "test_data"
DATA_REAL = Path(__file__).resolve().parent / "data_real"
DATA_GENOMES = Path(__file__).resolve().parent / "data_genomes"


@dataclass
class Workload:
    id: str
    description: str
    # Returns (target_path, query_path, meta_dict). Called lazily so synthetic
    # workloads only generate when actually run.
    resolve: Callable[[], tuple[Path, Path, dict]]
    tags: tuple[str, ...] = ()
    # lastz "actions" appended to each input as [a1,a2,...]. Most commonly
    # `multiple` (allow >1 sequence in the file). See `lastz --help=files`.
    target_actions: tuple[str, ...] = ()
    query_actions:  tuple[str, ...] = ()
    # Workload-specific lastz flags, applied before variant.extra_args. Use
    # this for things the workload itself requires (e.g. --format=maf when the
    # default lav format is incompatible with multi-sequence targets).
    lastz_args: tuple[str, ...] = ()


@dataclass
class Variant:
    id: str
    description: str
    # Extra args appended after the two input files. Should NOT include
    # --output= (the runner controls that).
    extra_args: list[str] = field(default_factory=list)


def _static(target_rel: str, query_rel: str) -> Callable[[], tuple[Path, Path, dict]]:
    def _resolve() -> tuple[Path, Path, dict]:
        t = LASTZ_TEST_DATA / target_rel
        q = LASTZ_TEST_DATA / query_rel
        if not t.exists() or not q.exists():
            raise FileNotFoundError(f"static workload missing inputs: {t} / {q}")
        meta = {"kind": "static", "target": str(t), "query": str(q)}
        return t, q, meta
    return _resolve


def _real(target_path: Path, query_path: Path) -> Callable[[], tuple[Path, Path, dict]]:
    def _resolve() -> tuple[Path, Path, dict]:
        if not target_path.exists() or not query_path.exists():
            raise FileNotFoundError(
                f"real-data workload missing inputs: {target_path} / {query_path}.\n"
                f"Run `python3 bench/fetch_real_data.py` to fetch them."
            )
        meta = {"kind": "real", "target": str(target_path), "query": str(query_path)}
        return target_path, query_path, meta
    return _resolve


def _synthetic(
    name: str,
    target_len: int,
    query_len: int | None = None,
    sub_rate: float = 0.05,
    indel_rate: float = 0.01,
    seed: int = 0xC0FFEE,
) -> Callable[[], tuple[Path, Path, dict]]:
    def _resolve() -> tuple[Path, Path, dict]:
        from gen_data import make_pair  # local import to keep import cost low
        out_dir = Path(__file__).parent / "data"
        t, q, meta = make_pair(
            out_dir=out_dir,
            base_name=name,
            target_len=target_len,
            query_len=query_len,
            sub_rate=sub_rate,
            indel_rate=indel_rate,
            seed=seed,
        )
        meta["kind"] = "synthetic"
        return t, q, meta
    return _resolve


# ---------------------------------------------------------------------------
# Workload registry
# ---------------------------------------------------------------------------

WORKLOADS: dict[str, Workload] = {}


def _register(w: Workload) -> None:
    if w.id in WORKLOADS:
        raise ValueError(f"duplicate workload id: {w.id}")
    WORKLOADS[w.id] = w


# Static workloads sourced from the lastz test_data directory.
_register(Workload(
    id="static.pseudocat_vs_pseudopig",
    description="Bundled lastz smoke test: pseudocat.fa vs pseudopig.fa (~19KB / ~69KB).",
    resolve=_static("pseudocat.fa", "pseudopig.fa"),
    tags=("smoke", "tiny"),
))

_register(Workload(
    id="static.pseudopig_self",
    description="pseudopig.fa vs pseudopig2.fa — near-self alignment, lots of HSPs.",
    resolve=_static("pseudopig.fa", "pseudopig2.fa"),
    tags=("smoke", "tiny"),
    # pseudopig.fa contains multiple sequences; lastz needs `multiple` to allow
    # it, and the default `lav` format rejects `multiple` so we emit `maf`.
    target_actions=("multiple",),
    lastz_args=("--format=maf",),
))

_register(Workload(
    id="static.reads_vs_pseudopig",
    description="sample_101s.fa (short reads) vs pseudopig.fa — read-mapping shape.",
    resolve=_static("pseudopig.fa", "sample_101s.fa"),
    tags=("smoke", "tiny", "reads"),
    target_actions=("multiple",),
    lastz_args=("--format=maf",),
))

# Real-data: chromosome-scale slices for SegAlign-comparable benchmarks.
# These point at ~60–130 Mbp single-chromosome FASTAs that bench/fetch_genomes.py
# extracts from the UCSC 2bit assemblies. Soft-masking is preserved.
_register(Workload(
    id="real.galGal6_chrZ_self",
    description=("Chicken (galGal6) chrZ vs itself, ~82 Mbp x 82 Mbp. Self-vs-self: "
                 "trivial diagonal plus any paralogous segments / unmasked repeats. "
                 "First real chromosome-scale calibration."),
    resolve=_real(
        target_path=DATA_GENOMES / "galGal6.chrZ.fa",
        query_path=DATA_GENOMES / "galGal6.chrZ.fa",
    ),
    tags=("real", "genome", "chromosome", "self", "bird"),
))

_register(Workload(
    id="real.galGal6_vs_taeGut2_chrZ",
    description=("Chicken (galGal6) chrZ vs Zebra Finch (taeGut2) chrZ, "
                 "~82 Mbp x ~72 Mbp. Bird-vs-bird, ~100 Mya divergence — "
                 "expect strong syntenic block alignments."),
    resolve=_real(
        target_path=DATA_GENOMES / "galGal6.chrZ.fa",
        query_path=DATA_GENOMES / "taeGut2.chrZ.fa",
    ),
    tags=("real", "genome", "chromosome", "cross", "bird"),
))

_register(Workload(
    id="real.galGal6_vs_taeGut2_chrZ_10mb",
    description=("First 10 Mbp of chicken chrZ vs first 10 Mbp of zebra "
                 "finch chrZ. Iteration-friendly subset of the bird-Z "
                 "cross workload (~1 min wall) for profiler development. "
                 "Same shape as the full chrZ pair: real DNA, real "
                 "soft-masking, real cross-species divergence — just "
                 "1/(8 x 7) of the seed-search space."),
    resolve=_real(
        target_path=DATA_GENOMES / "galGal6.chrZ_0_10mb.fa",
        query_path=DATA_GENOMES / "taeGut2.chrZ_0_10mb.fa",
    ),
    tags=("real", "genome", "subset", "cross", "bird", "profile"),
))

_register(Workload(
    id="real.hg38_vs_mm10",
    description=("Human (hg38) chr19 vs Mouse (mm10) chr10, ~58 Mbp x ~130 Mbp. "
                 "Mostly-syntenic mammal pair, ~90 Mya divergence. Heavier than "
                 "the bird pair due to higher repeat content."),
    resolve=_real(
        target_path=DATA_GENOMES / "hg38.chr19.fa",
        query_path=DATA_GENOMES / "mm10.chr10.fa",
    ),
    tags=("real", "genome", "chromosome", "cross", "mammal"),
))

# Real-data: cross-species β-actin mRNA matrix (human / chimp / fruit fly).
# Each cell is target=<species_t> × query=<species_q>. Sequences are RefSeq
# mRNAs of comparable size (~1.8–1.9 kb). See bench/fetch_real_data.py.
ACTB_SPECIES = ("human", "chimp", "fly")
for _t in ACTB_SPECIES:
    for _q in ACTB_SPECIES:
        _register(Workload(
            id=f"real.actb_{_t}_vs_{_q}",
            description=(f"β-actin mRNA: {_t} (target) vs {_q} (query). "
                         f"Tiny (~1.8 kb). Used to study how alignment shape "
                         f"changes with evolutionary distance."),
            resolve=_real(
                target_path=DATA_REAL / f"actb_{_t}.fa",
                query_path=DATA_REAL / f"actb_{_q}.fa",
            ),
            tags=("real", "actb", "tiny", "matrix"),
        ))


# Synthetic scaling sweep. Sizes chosen to span ~3 orders of magnitude with
# reasonable wall-clock budget on a single core.
for size_kb, tag in [
    (10, "tiny"),
    (100, "small"),
    (1_000, "medium"),
    (10_000, "large"),
]:
    _register(Workload(
        id=f"synth.uniform_{size_kb}kb",
        description=f"Synthetic random DNA, {size_kb} kbp target vs mutated query "
                    f"(sub=5%, indel=1%).",
        resolve=_synthetic(
            name=f"uniform_{size_kb}kb",
            target_len=size_kb * 1_000,
        ),
        tags=("synth", "scaling", tag),
    ))


# ---------------------------------------------------------------------------
# Variant registry
# ---------------------------------------------------------------------------

VARIANTS: dict[str, Variant] = {
    "default": Variant(
        id="default",
        description="lastz with all defaults (seed-and-extend + gapped, no chaining).",
        extra_args=[],
    ),
    # SegAlign uses lastz's defaults (seed=12of19, step=1, xdrop=910,
    # hspthresh=3000, ydrop≈9430, gappedthresh=hspthresh) — confirmed against
    # gsneha26/SegAlign src/main.cpp. The only practical difference vs our
    # `default` is the output format: SegAlign writes `maf-` (MAF without the
    # `# lastz ...` header). Keep this as the canonical SegAlign-comparable
    # variant so the runner output is byte-for-byte comparable.
    "segalign_default": Variant(
        id="segalign_default",
        description="Mirror of SegAlign's default invocation — lastz defaults + --format=maf-.",
        extra_args=["--format=maf-"],
    ),
    "nogapped": Variant(
        id="nogapped",
        description="Skip gapped extension (isolates seeding + ungapped extension).",
        extra_args=["--nogapped"],
    ),
    "chain": Variant(
        id="chain",
        description="Default + chaining stage on top.",
        extra_args=["--chain"],
    ),
    "step20": Variant(
        id="step20",
        description="--step=20 — sparser seeds, faster but lower sensitivity.",
        extra_args=["--step=20"],
    ),
}


# ---------------------------------------------------------------------------
# Suites
# ---------------------------------------------------------------------------

SUITES: dict[str, list[tuple[str, str]]] = {
    # Quick sanity sweep — runs in seconds, exercises every variant once.
    "smoke": [
        ("static.pseudocat_vs_pseudopig", "default"),
        ("static.pseudocat_vs_pseudopig", "nogapped"),
        ("static.pseudocat_vs_pseudopig", "chain"),
        ("static.pseudopig_self",        "default"),
        ("static.reads_vs_pseudopig",    "default"),
    ],
    # Synthetic scaling along input size, default flags only.
    "scaling": [
        ("synth.uniform_10kb",     "default"),
        ("synth.uniform_100kb",    "default"),
        ("synth.uniform_1000kb",   "default"),
        ("synth.uniform_10000kb",  "default"),
    ],
    # Same scaling but with gapped extension off — diff vs scaling tells us
    # how much of the runtime that stage owns.
    "scaling_nogapped": [
        ("synth.uniform_10kb",     "nogapped"),
        ("synth.uniform_100kb",    "nogapped"),
        ("synth.uniform_1000kb",   "nogapped"),
        ("synth.uniform_10000kb",  "nogapped"),
    ],
    # 3x3 cross-species β-actin mRNA matrix at default sensitivity. Expect
    # human↔chimp to give a single ~full-length alignment, and human/chimp↔fly
    # to give few or no HSPs (lastz defaults are tuned for ~70%+ identity).
    "actb_matrix": [
        (f"real.actb_{t}_vs_{q}", "default")
        for t in ACTB_SPECIES for q in ACTB_SPECIES
    ],
    # First real chromosome-scale calibration: bird Z chromosome self-vs-self.
    # Mirrors a SegAlign-style invocation. Single rep, expected ~5–30 min.
    "bird_z_self": [
        ("real.galGal6_chrZ_self", "segalign_default"),
    ],
    # Same workload, but break out by stage. `nogapped` isolates seed+ungapped
    # so subtracting it from `default` gives the gapped-extension cost.
    "bird_z_self_breakdown": [
        ("real.galGal6_chrZ_self", "default"),
        ("real.galGal6_chrZ_self", "nogapped"),
    ],
    # Bird-vs-bird cross-species: chicken Z vs zebra finch Z (~100 Mya
    # divergence). Z is sex-chromosome and well-conserved in birds, so we
    # expect strong syntenic blocks. This is the SegAlign-comparable shape:
    # two unrelated genomes, real masking, real evolutionary signal.
    "bird_z_cross": [
        ("real.galGal6_vs_taeGut2_chrZ", "segalign_default"),
    ],
    # Cross + breakdown: same pair under default and nogapped to back out the
    # gapped-extension share for a realistic (not self-trivial) workload.
    "bird_z_cross_breakdown": [
        ("real.galGal6_vs_taeGut2_chrZ", "default"),
        ("real.galGal6_vs_taeGut2_chrZ", "nogapped"),
    ],
    # Iteration-friendly subset (10 Mbp x 10 Mbp). For profiler development,
    # eu-stack sampling, and rdtsc substage breakdown — fast enough to iterate
    # but the same shape as the full chr workload.
    "bird_z_cross_10mb": [
        ("real.galGal6_vs_taeGut2_chrZ_10mb", "segalign_default"),
    ],
}

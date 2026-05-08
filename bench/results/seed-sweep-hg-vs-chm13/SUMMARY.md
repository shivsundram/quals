Within-species seed sweep on hg38 vs T2T-CHM13 chr1:50M-60M
============================================================

Workload: 10 Mbp × 10 Mbp at chr1:50,000,000-60,000,000.
Both assemblies are human; expected pairwise identity in single-copy
regions ~99.5-99.9% (CHM13 is a complete hydatidiform mole assembly,
hg38 is a mosaic; differences are real polymorphisms + repeat
instance variation, not species divergence).

Single-thread, pinned to CPU 0, lastz_TS (rdtsc-instrumented).

Headline metrics
----------------

variant                              wall    seed_hit  gapped_ext  HSPs   total_bp_aligned   mean_HSP_len
cross_species_default (12of19)       47.6 s    18.2 s    28.5 s    2474       12.44 Mbp        5,029 bp
bigger_seed_match15  (match15 step=1) 26.8 s     0.83 s   24.9 s    1798       12.12 Mbp        6,743 bp
hg_recommended (match15 step=20      21.5 s     0.37 s   20.5 s     943       11.33 Mbp       12,014 bp
              + notransition)

Wall speedup: bigger_seed = 1.78x, hg_recommended = 2.21x vs default.

Stage-shift: who's the bottleneck?
-----------------------------------

12of19 (cross-species default):       seed_hit = 38%, gapped = 60%  (SHARED)
match15:                              seed_hit =  3%, gapped = 93%  (gapped dominant)
match15+step20:                       seed_hit =  2%, gapped = 95%  (gapped dominant)

This is a major regime change vs. the bird-Z 10 Mbp cross-species run,
where seed_hit_search was 88% of total and gapped_extension was 10%.
At high identity with appropriate seed tuning, the bottleneck FLIPS:
gapped extension becomes 90%+ of runtime. SegAlign's GPU port of both
seed_hit_search AND gapped extension is necessary for human-vs-human
work; just doing seed_hit_search would leave 95% of the runtime on the
CPU.

Sensitivity check: HSP merging vs. true homology loss
-----------------------------------------------------

The HSP count drops 62% (2474 -> 943) but total-bp-covered drops only
9% (12.44 -> 11.33 Mbp). Mean HSP length more than doubles (5 KB ->
12 KB). Maximum HSP length is identical across all three (~448 KB).

So "fewer HSPs" is misleading -- the missing HSPs aren't missing
homology, they're short fragments that get reported as single longer
HSPs when the seed grid is sparser. Total-bp-aligned is the right
sensitivity metric, and at 9% loss for ~50% wall savings this is a
strong Pareto improvement.

The 9% bp loss is itself probably mostly real homologs at the lowest
identity tail (e.g. partially-diverged repeat instances) where 15
contiguous matching positions never exist.

Why does seed_hit_search drop ~50x but wall time only drops ~2x?
-----------------------------------------------------------------

Because for cross-species defaults on this workload, seed_hit_search
was already ONLY 38% of total. The other 60% is gapped extension,
which is set by HSP count and HSP length, not by seed-table machinery.
Going from 12of19 -> match15+step20 drops gapped from 28.5s -> 20.5s
(28% reduction) -- because there are fewer surviving HSP candidates,
not because each HSP is faster to process.

To see the seed-machinery wins translate to bigger wall speedups, you
either need:
  (a) workloads where seed_hit_search is a larger share of total
      (cross-species, very deep divergence -- bird-Z 10 Mbp had
      seed_hit_search at 88%), OR
  (b) GPU acceleration of gapped extension too, so the new bottleneck
      gets attacked.

rdtsc instrumentation note
--------------------------

The rdtsc counters in find_table_matches show 0 cycles for the match15
runs, because match15 has bit-weight 30 (> 24-bit max index), so lastz
takes the "overweight seed" path through find_table_matches_resolve()
which is a separate function we haven't instrumented. The xdrop and
reporter timers still fire correctly. To get full rdtsc coverage on
match15-class seeds we'd need to add the same instrumentation block
to find_table_matches_resolve. Wall-clock data is unaffected.

Also: the dedup short-circuit fires 96% of pfs calls on match15 vs 1%
on 12of19. That's because contiguous seeds funnel repeats into a
single bucket (no don't-care positions to scatter them across), so
the diagEnd[hDiag] check is doing real work again -- exactly the
opposite of the bird-Z observation. Confirms that dedup's value
depends heavily on seed shape, not just seed weight.

Conclusion
----------

For human-vs-human at vertebrate-genome scale:
  - seed_hit_search bottleneck (the canonical "lastz is slow because
    of seed lookup" claim) does not hold; bigger seeds + step>1 cut
    it to under 2% of runtime.
  - The ACTUAL bottleneck becomes gapped extension. SegAlign's GPU
    port of gapped extension is what matters here.
  - Sensitivity loss from bigger/sparser seeds is small (~9% bp loss
    for ~55% wall savings) -- well worth it for any practical
    pangenome / assembly-vs-assembly workflow.
  - Pure "bigger seed" wins are real but capped at ~2x wall speedup
    on this workload because of the gapped-extension floor. To go
    further, gapped extension must also be accelerated.

This is the regime change SegAlign was designed for vertebrate-scale
cross-species (bird-Z), where seed_hit_search is 88% of total. In
within-species mode, their seed_hit_search wins matter less; their
gapped extension wins matter more.

Seed-weight sweep on bird-Z 10 Mbp × 10 Mbp (galGal6.chrZ vs taeGut2.chrZ)
=========================================================================

Single-thread run on the head node, lastz_TS (rdtsc-instrumented), pinned
to CPU 0. All other parameters identical (segalign-style: HOXD, ydrop=15000,
gappedthresh=3000, hspthresh=3000, --format=maf-).

Headline metrics
----------------

variant         wall    total   seed_hit   chain_iters   xdrop_calls   xdrop_no_score   HSPs   sens vs 12of19
12of19 (def)    62.3 s  62.3 s   55.2 s      164.8 M       163.0 M       99.999 %       510    100.0 %
14of22 (heavy)  17.1 s  17.1 s    9.9 s       13.2 M        13.2 M        99.988 %       490     96.1 %  (lost 20)
match12 (cont)  10.4 s  10.4 s    4.3 s       18.7 M        13.6 M        99.990 %       436     85.5 %  (lost 74)

Speedups: 14of22 = 3.6×,  match12 = 6.0×.

rdtsc substage breakdown (per-call cycle averages)
--------------------------------------------------

substage              12of19       14of22       match12     interpretation
ftm.chain_walk        125.5        49.7         125.8       per-iter pointer-chase; same for 12of19/match12 (same table size); much
                                                            cheaper for 14of22 because most lookups hit empty buckets and never enter
                                                            the loop, so only short hot chains contribute.
ftm.processor_cb      495.6        554.6        288.4       wraps pfs+xdrop+reporter; tracks xdrop cost mostly.
pfs.dedup_skip cnt     1.85 M       0.07 M       5.08 M     27% of pfs calls in match12 hit dedup vs 1% in 12of19!
pfs.dedup_skip frac    1.1 %        0.5 %       27.2 %      contig seed exposes repeats that spaced seed scatters across buckets.
xdrop_extend cyc/call  424.5        497.8        277.1       contig hits die faster (low-conservation random matches).
reporter_calls         1665         1588         1298        ungapped HSPs that survived xdrop (vs MAF blocks below).

Final HSPs (MAF blocks)
-----------------------

12of19: 510 blocks (1.44 MB MAF output)
14of22: 490 blocks (1.41 MB)        — lost 20 (3.9%) from heavier specificity
match12: 436 blocks (1.38 MB)       — lost 74 (14.5%) from no don't-care positions

Key takeaways
-------------

1. 14of22 is a near-Pareto win for cross-species at ~75% identity:
   3.6× speedup for 4% sensitivity loss. SegAlign keeps 12of19 as the safe
   default for distant species, but 14of22 is the right choice when species
   are close enough that you can afford slightly less sensitivity.

2. Contiguous match12 loses 14.5% of HSPs — that's real homology missed.
   The 6× speedup is not worth it for a quals-grade aligner. Confirms that
   spaced seeds are doing real work: at the same weight (12 required positions)
   they're 10% more sensitive than contiguous.

3. Diagonal dedup is dramatically more useful for contig seeds than for spaced
   seeds. match12 hits 27% dedup short-circuits vs only 1.1% for 12of19. The
   spaced seed's don't-care positions scatter repeat hits across different
   buckets, doing most of dedup's job before lookup. SegAlign throwing dedup
   away on the GPU costs even less than we measured: 0.5% of pfs calls on
   14of22, 1.1% on 12of19. Biggest cost would be on contig seeds (27%) which
   nobody actually uses.

4. The chain-walk DRAM-latency story (125 cyc/iter) is INDEPENDENT of seed
   weight: same on 12of19 and match12 (both use 4^12 buckets). Going to 14of22
   shrinks the chain-walk per-iter cost only because most lookups hit empty
   buckets and short-circuit before entering the loop — not because the
   underlying memory pattern changed. The dense-array layout fix would still
   apply at any seed weight.

5. The xdrop noScore rate stays at 99.99%+ for ALL three variants. Heavier
   seeds reduce the absolute number of failed extensions (163M → 13M), but
   the success rate is fundamentally seed-independent: ~1% of x-drop calls
   ever produce a real HSP. The GPU's "embarrassingly parallel failed
   extension" win applies at any seed weight.

Conclusion: seed weight is a sensitivity/speed dial, NOT a fix for the
algorithmic problems we identified. The real wins remain (a) dense-array
position table for chain-walk DRAM-latency, (b) GPU x-drop for the failed-
extension parallelism. Both compound with whichever seed weight you pick.

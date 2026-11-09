#!/usr/bin/env python3
"""Shared unmappable intervals: regions with no recombination-map data in ANY
population.

Why not a per-population gap threshold: map density varies a lot between
populations (chr2 has 136,733 map rows for YRI but 88,046 for CEU), so
"interval wider than X" measures the MAP, not the GENOME. Applied per
population it dropped 0.003% of YRI sites and 0.253% of CHB sites, which
would leave the four core populations on different site lists -- and would
discard real, mappable regions in the sparser maps, where linear
interpolation across a wide interval is the ordinary assumption of constant
rate and is exactly what every genetic map does.

A centromere is an assembly feature and is identical in every population, so
the exclusion must be too. An interval is excluded only if it is a gap in
EVERY population's map, i.e. the intersection of the per-population gaps.
That privileges no population and leaves all four with one site list.

Writes gaps.tsv: chrom, start, end, n_pop (always all), width.
"""
import argparse
import os
import sys

MAP_FMT = "{d}/{pop}/{pop}_map_chr_{chrom}.txt"


def gaps_for(path, min_gap):
    """-> sorted [(start, end)] of consecutive map rows further than min_gap."""
    pos = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 3 and p[0].isdigit():
                pos.append(int(p[0]))
    pos.sort()
    out = [(pos[i], pos[i + 1]) for i in range(len(pos) - 1)
           if pos[i + 1] - pos[i] > min_gap]
    # outside the mapped range is unmappable too, and is handled separately
    # per population in fix_cm.py (a few hundred sites genome-wide).
    return out, (pos[0], pos[-1])


def intersect(intervals):
    """Intersection of several interval lists (each sorted, disjoint)."""
    cur = intervals[0]
    for nxt in intervals[1:]:
        out, i, j = [], 0, 0
        while i < len(cur) and j < len(nxt):
            lo = max(cur[i][0], nxt[j][0])
            hi = min(cur[i][1], nxt[j][1])
            if lo < hi:
                out.append((lo, hi))
            if cur[i][1] < nxt[j][1]:
                i += 1
            else:
                j += 1
        cur = out
    return cur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pops", default="YRI,CEU,CHB,BEB")
    ap.add_argument("--map-dir", required=True, help="<dir>/<POP>/<POP>_map_chr_<N>.txt (convert_pyrho_maps.py)")
    ap.add_argument("--min-gap", type=int, default=1_000_000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pops = a.pops.split(",")

    rows = []
    for chrom in range(1, 23):
        per_pop, ranges = [], []
        for pop in pops:
            f = MAP_FMT.format(d=a.map_dir, pop=pop, chrom=chrom)
            if not os.path.exists(f):
                sys.exit(f"ERROR: map not found: {f}")
            g, rng = gaps_for(f, a.min_gap)
            per_pop.append(g)
            ranges.append(rng)
        shared = intersect(per_pop)
        for lo, hi in shared:
            rows.append((chrom, lo, hi, len(pops), hi - lo))
        print(f"chr{chrom}: per-pop gaps {[len(g) for g in per_pop]} -> "
              f"shared {len(shared)} covering {sum(h - l for l, h in shared) / 1e6:.1f} Mb",
              flush=True)

    with open(a.out, "w") as fh:
        fh.write("chrom\tstart\tend\tn_pop\twidth\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"\n{len(rows)} shared unmappable intervals -> {a.out}")
    print(f"total {sum(r[4] for r in rows) / 1e6:.1f} Mb excluded")


if __name__ == "__main__":
    main()

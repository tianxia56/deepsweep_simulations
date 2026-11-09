#!/usr/bin/env python3
"""Genome-wide exact-k reference tables for one target population.

Pools the per-chromosome (k, value) arrays written by emp_stats.py over all 22
autosomes and writes one table per statistic in the SAME format and with the
same rules as the simulated bins (cms_stats.build_bins): one row per exact
derived count k, columns k, n_hap, mean, std, count, and a bin holding fewer
than MIN_SITES_PER_BIN sites written NaN rather than given an unreliable sigma.

Two things differ from the simulated reference, both by design:
  * the pool is GENOME-WIDE EMPIRICAL, not neutral, so it contains whatever
    real selection the genome carries and is shifted by the signal being
    measured. This is the user's specified setup; it is recorded in the
    MANIFEST so no later reader mistakes it for a neutral null.
  * k runs over the POPULATION'S OWN n_hap (YRI 216, CEU 198, CHB 206,
    BEB 172). Exact-k tables are only valid at the sample size they were
    built at, which is why no table is shared between populations.

Writes bins/<target>/{nsl,ihh12,ihs_binom,delihh_binom}_bin.csv,
xpehh_<target>_vs_<alt>_bin.csv and MANIFEST.txt.
"""
import argparse
import glob
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "pipeline"))   # cms_stats
import cms_stats as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chroms", default=",".join(str(c) for c in range(1, 23)))
    a = ap.parse_args()
    chroms = [int(c) for c in a.chroms.split(",")]

    files = []
    for c in chroms:
        f = f"{a.raw_dir}/{a.target}.{c}.kvals.npz"
        if not os.path.exists(f):
            sys.exit(f"ERROR: missing {f} - run emp_stats.py for every chromosome")
        files.append((c, f))

    pools, n_hap = {}, None
    for c, f in files:
        d = np.load(f)
        nh = int(d["n_hap"])
        if n_hap is None:
            n_hap = nh
        elif nh != n_hap:
            sys.exit(f"ERROR: {f} has n_hap {nh}, expected {n_hap}")
        for key in d.files:
            if not key.endswith("__k"):
                continue
            stat = key[:-3]
            pools.setdefault(stat, [[], []])
            pools[stat][0].append(d[key])
            pools[stat][1].append(d[f"{stat}__v"])

    os.makedirs(a.out_dir, exist_ok=True)
    rows = []
    for stat in sorted(pools):
        k = np.concatenate(pools[stat][0])
        v = np.concatenate(pools[stat][1])
        # build_bins takes DAF; k/n_hap is exact and round-trips through
        # count_from_daf, so the table is keyed on the same k.
        tab = C.build_bins(k / n_hap, v, n_hap,
                           out_csv=f"{a.out_dir}/{stat}_bin.csv")
        usable = int(tab["mean"].notna().sum())
        rows.append((stat, len(v), len(tab), usable, int(tab["count"].min()),
                     int(tab["count"].max())))
        print(f"{stat:<28} {len(v):>12,} sites  {usable}/{len(tab)} usable bins",
              flush=True)

    with open(f"{a.out_dir}/MANIFEST.txt", "w") as fh:
        fh.write(f"built          : {datetime.now():%F %T}\n")
        fh.write(f"target         : {a.target}\n")
        fh.write(f"n_hap          : {n_hap}\n")
        fh.write(f"chromosomes    : {len(chroms)} autosomes\n")
        fh.write("reference      : GENOME-WIDE EMPIRICAL, *not* neutral. Contains\n")
        fh.write("                 whatever real selection the genome carries, so it\n")
        fh.write("                 is shifted by the signal being measured. It is NOT\n")
        fh.write("                 equivalent to the simulated neutral null.\n")
        fh.write("site list      : target + its three outgroups -- a position is\n")
        fh.write("                 kept when at least one of those four is\n")
        fh.write("                 polymorphic there, intersected with the positions\n")
        fh.write("                 all four carry; applied AFTER the tools ran\n")
        fh.write(f"binning        : one bin per exact derived count k, 0..{n_hap}\n")
        fh.write(f"MIN_SITES_PER_BIN : {C.MIN_SITES_PER_BIN} (thinner bins -> NaN)\n")
        fh.write(f"selscan        : {' '.join(C.SELSCAN_FLAGS)}\n")
        fh.write(f"ihsbin         : {' '.join(C.IHSBIN_FLAGS)} --binom (delihh SIGNED)\n")
        fh.write(f"xpehhbin       : {' '.join(C.XPEHHBIN_FLAGS)}\n")
        fh.write(f"XPEHH_KEEP_TARGET_MONOMORPHIC : {C.XPEHH_KEEP_TARGET_MONOMORPHIC}\n")
        fh.write("FST            : Weir-Hill; fixed-same pairs NaN, excluded from mean\n")
        fh.write("cM             : linear interpolation (annotate_tped.py); shared unmappable\n")
        fh.write("                 intervals dropped (build_gap_mask.py gaps.tsv)\n")
        for t in C.tool_versions().items():
            fh.write(f"tool           : {t[0]}: {t[1]}\n")
        fh.write("\nstat\tn_sites\tbins\tusable\tmin_count\tmax_count\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"\n-> {a.out_dir}/MANIFEST.txt")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Normalize one (target, chromosome) against the genome-wide exact-k tables
and write the final per-chromosome table.

Uses cms_stats.normalize, the same function the simulated corpus uses, so the
z-scores are produced by identical code. max_xpehh is the row-wise max of the
per-pair NORMALIZED scores, joined on variant ID, exactly as 02 does it.

Output columns are the simulated schema minus the sweep columns (which have no
empirical meaning), plus chr:
  chr pos ID mean_fst deldaf daf maf norm_nsl norm_ihh12 norm_ihs_binom
  norm_delihh_binom max_xpehh
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "pipeline"))   # cms_stats
sys.path.insert(0, HERE)
import cms_stats as C  # noqa: E402
from emp_stats import outgroups_for  # noqa: E402  single definition

CORE4 = ["YRI", "CEU", "CHB", "BEB"]
OUT_COLS = ["chr", "pos", "ID", "mean_fst", "deldaf", "daf", "maf",
            "norm_nsl", "norm_ihh12", "norm_ihs_binom", "norm_delihh_binom",
            "max_xpehh"]


def log(m):
    print(f"[{time.strftime('%F %T')}] {m}", flush=True)


def read_manifest_n_hap(bins_dir):
    for line in open(f"{bins_dir}/MANIFEST.txt"):
        if line.startswith("n_hap"):
            return int(line.split(":", 1)[1])
    sys.exit(f"ERROR: {bins_dir}/MANIFEST.txt records no n_hap")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--chrom", type=int, required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--bins-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--state-dir", required=True)
    ap.add_argument("--outgroups", default=None,
                    help="comma list; default = the core populations outside "
                         "the target's superpopulation (emp_stats.outgroups_for)")
    a = ap.parse_args()

    outg = (a.outgroups.split(",") if a.outgroups
            else outgroups_for(a.target))
    if len(outg) != 3 or a.target in outg:
        sys.exit(f"ERROR: bad outgroup set for {a.target}: {outg}")
    state = os.path.join(a.state_dir, f"{a.target}__{a.chrom}.json")
    if os.path.exists(state):
        log(f"{a.target} chr{a.chrom}: already done")
        return
    for d in (a.out_dir, a.state_dir):
        os.makedirs(d, exist_ok=True)

    n_hap = read_manifest_n_hap(a.bins_dir)
    raw = f"{a.raw_dir}/{a.target}.{a.chrom}.raw.tsv"
    df = pd.read_csv(raw, sep="\t", dtype={"ID": str})
    kv = np.load(f"{a.raw_dir}/{a.target}.{a.chrom}.kvals.npz")
    if int(kv["n_hap"]) != n_hap:
        sys.exit(f"ERROR: raw n_hap {int(kv['n_hap'])} != bins n_hap {n_hap}")
    t0 = time.time()

    # ---- one-population statistics ---------------------------------------
    for src, dst in (("nsl", "norm_nsl"), ("ihh12", "norm_ihh12"),
                     ("ihs_binom", "norm_ihs_binom"),
                     ("delihh_binom", "norm_delihh_binom")):
        tab = pd.read_csv(f"{a.bins_dir}/{src}_bin.csv")
        df[dst] = C.normalize(df[src].to_numpy(dtype=float),
                              df["daf"].to_numpy(), tab, n_hap)

    # ---- XP-EHH: normalize each pair, then row-wise max on variant ID ----
    # The raw per-pair outputs are not kept (they are 10x the final table), so
    # the pair values are recovered from the kvals arrays, which were written
    # already restricted to the kept sites and in that order.
    xp_cols = []
    for i, alt in enumerate(outg):
        pair = f"{a.target}_vs_{alt}"
        f = f"{a.raw_dir}/{a.target}.{a.chrom}.xp_{pair}.tsv"
        if not os.path.exists(f):
            sys.exit(f"ERROR: missing per-pair XP-EHH table {f}")
        d = pd.read_csv(f, sep="\t", dtype={"ID": str})
        tab = pd.read_csv(f"{a.bins_dir}/xpehh_{pair}_bin.csv")
        col = f"_xp{i}"
        d[col] = C.normalize(d["xpehh"].to_numpy(dtype=float),
                             d["daf"].to_numpy(), tab, n_hap)
        df = df.merge(d[["ID", col]], on="ID", how="left")
        xp_cols.append(col)
    df["max_xpehh"] = df[xp_cols].max(axis=1) if xp_cols else np.nan
    df = df.drop(columns=xp_cols)

    df["chr"] = a.chrom
    out = df.reindex(columns=OUT_COLS).round(6)
    tmp = f"{a.out_dir}/{a.target}.{a.chrom}.tsv.tmp"
    out.to_csv(tmp, sep="\t", index=False, na_rep="NA")
    os.replace(tmp, f"{a.out_dir}/{a.target}.{a.chrom}.tsv")

    import json
    info = dict(target=a.target, chrom=a.chrom, n_hap=n_hap, rows=len(out),
                outgroups=outg, secs=round(time.time() - t0, 1),
                na={c: int(out[c].isna().sum()) for c in OUT_COLS[3:]},
                z={c: dict(mean=round(float(out[c].mean(skipna=True)), 4),
                           sd=round(float(out[c].std(skipna=True)), 4))
                   for c in OUT_COLS[7:]})
    with open(state + ".tmp", "w") as fh:
        json.dump(info, fh, indent=1)
    os.replace(state + ".tmp", state)
    log(f"{a.target} chr{a.chrom}: {len(out):,} rows, {info['secs']} s")


if __name__ == "__main__":
    main()

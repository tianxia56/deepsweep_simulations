#!/usr/bin/env python3
"""Raw component statistics for one (target population, chromosome).

Every tool, flag and formula is the simulation pipeline's, taken from
cms_stats so the two cannot drift: selscan v3.0.0 --maf 0.05 (nSL, iHH12),
ihsbin --binom --minmaf 0.05 (iHS and SIGNED dihh), xpehhbin --minmaf 0,
Weir-Hill FST with fixed-same pairs NaN and excluded from mean_fst,
deldaf = daf_target - mean(daf_outgroups), derived = 1.

Two deliberate departures, both specified by the user:
  * the null is GENOME-WIDE EMPIRICAL, built per population at its OWN n_hap
    (YRI 216, CEU 198, CHB 206, BEB 172), not the neutral simulations. This
    stage only writes the raw values and their exact derived count k;
    build_bins.py pools them genome-wide and normalize.py applies them.
  * the site list is the TARGET AND ITS THREE OUTGROUPS: a position is kept
    when at least one of those four is polymorphic there, intersected with the
    positions all four actually carry, which FST and deldaf need. Applied
    AFTER every tool has run, never before: thinning the grid would change the
    EHH integral.

The one-population statistics run on the target's OWN full site list, which
is the analogue of the simulated target TPED (every site in the genealogy,
target-monomorphic ones included). XP-EHH runs on each pair's harmonised
intersection, because xpehhbin needs both populations on one grid.

Writes out/<target>.<chrom>.raw.tsv (one row per kept site),
out/<target>.<chrom>.kvals.npz (exact k + value per statistic, for binning)
and out/<target>.<chrom>.xp_<pair>.tsv (ID, daf, xpehh per pair, kept sites
only). The per-pair tables are what normalize.py needs: max_xpehh is the max
of the per-pair NORMALIZED scores, so each pair has to be z-scored against
its own table before the max, and that merge is on variant ID.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "pipeline"))   # cms_stats
import cms_stats as C  # noqa: E402

CORE4 = ["YRI", "CEU", "CHB", "BEB"]
HARMONIZE = os.path.join(HERE, "harmonize.py")

# Outgroups = the core populations NOT in the target's superpopulation, so
# every target gets exactly three (calc_fst_daf_deldaf.py:107). That keeps
# max_xpehh a max over three normalised pairs for all 22 populations, as it
# is for core4 and for the 4-population simulated models. The SITE LIST is
# these same four populations: a position is kept when at least one of them is
# polymorphic there.
SUPERPOP = {
    "ACB": "AFR", "ASW": "AFR", "ESN": "AFR", "GWD": "AFR", "LWK": "AFR",
    "MSL": "AFR", "YRI": "AFR",
    "CEU": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR", "TSI": "EUR",
    "CDX": "EAS", "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "KHV": "EAS",
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
}
CORE_BY_SUPERPOP = {"AFR": "YRI", "EUR": "CEU", "EAS": "CHB", "SAS": "BEB"}


def outgroups_for(target):
    sp = SUPERPOP[target]
    return [c for s, c in CORE_BY_SUPERPOP.items() if s != sp]


def log(m):
    print(f"[{time.strftime('%F %T')}] {m}", flush=True)


def common_positions(parsed):
    """Intersect the populations' position arrays; -> {pop: index array}."""
    common = None
    for p in parsed:
        common = parsed[p][1] if common is None else np.intersect1d(common, parsed[p][1])
    idx = {}
    for p in parsed:
        pos = parsed[p][1]
        order = np.argsort(pos, kind="stable")
        hit = order[np.searchsorted(pos, common, sorter=order)]
        if not np.array_equal(pos[hit], common):
            sys.exit(f"ERROR: position intersection failed for pop {p}")
        idx[p] = hit
    return common, idx


def run_xpehh(target, alt, chrom, hmap_dir, tmp, out_prefix):
    """Harmonise the pair onto its shared grid, then xpehhbin. -> output path."""
    base = f"{tmp}/{target}_vs_{alt}.{chrom}"
    args = [f"{hmap_dir}/{target}.{chrom}.map", f"{hmap_dir}/{target}.{chrom}.hap",
            f"{hmap_dir}/{alt}.{chrom}.map", f"{hmap_dir}/{alt}.{chrom}.hap",
            f"{base}.map", f"{base}.hapA", f"{base}.hapB", f"{base}.sites.npz"]
    subprocess.run([sys.executable, HARMONIZE, *args], check=True)
    out = f"{out_prefix}_{target}_vs_{alt}.xpehh.out"
    C.run_xpehhbin(f"{base}.hapA", f"{base}.hapB", f"{base}.map", out)
    for f in (f"{base}.map", f"{base}.hapA", f"{base}.hapB", f"{base}.sites.npz"):
        os.remove(f)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(SUPERPOP))
    ap.add_argument("--chrom", type=int, required=True)
    ap.add_argument("--tped-dir", required=True)
    ap.add_argument("--hmap-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--state-dir", required=True)
    ap.add_argument("--keep-raw-xpehh", action="store_true")
    a = ap.parse_args()

    state = os.path.join(a.state_dir, f"{a.target}__{a.chrom}.json")
    if os.path.exists(state):
        log(f"{a.target} chr{a.chrom}: already done")
        return
    C.check_tools(["selscan", "ihsbin", "xpehhbin"])
    flavours = C.ihs_flavours(False)
    assert list(flavours) == ["binom"], flavours
    for d in (a.out_dir, a.state_dir):
        os.makedirs(d, exist_ok=True)
    alt_pops = outgroups_for(a.target)
    pops = [a.target] + alt_pops
    tmp = os.path.join(os.environ.get("TMPDIR", "/tmp"),
                       f"emp_{a.target}_{a.chrom}_{os.getpid()}")
    os.makedirs(tmp, exist_ok=True)
    t0 = time.time()

    try:
        # ---- parse all four populations (frequencies only) -------------------
        parsed = {}
        for p in pops:
            f = f"{a.tped_dir}/{p}.{a.chrom}.tped"
            if not os.path.exists(f):
                sys.exit(f"ERROR: missing {f}")
            parsed[p] = C.parse_tped(f)
            log(f"  {p}: {len(parsed[p][1]):,} sites, n_hap {parsed[p][3]}")
        n_hap = parsed[a.target][3]
        # The haplotype statistics are joined on ID, so a repeated ID would
        # make a Cartesian product and attach values from the wrong site.
        for p in pops:
            ids = parsed[p][0]
            if len(np.unique(ids)) != len(ids):
                sys.exit(f"ERROR: {p} chr{a.chrom} has duplicate variant IDs; "
                         f"rerun annotate_tped.py (it renames \".\" and asserts uniqueness)")

        common, idx = common_positions(parsed)
        log(f"  four-population common positions: {len(common):,}")
        sub = {p: (parsed[p][0][idx[p]], parsed[p][1][idx[p]],
                   parsed[p][2][idx[p]], parsed[p][3]) for p in pops}
        for p in pops:
            if p != a.target and not np.array_equal(sub[p][0], sub[a.target][0]):
                sys.exit(f"ERROR: rsID mismatch between {a.target} and {p} at "
                         f"shared positions")
        table = C.sfs_table(sub, a.target, alt_pops)

        # ---- one-population statistics, on the target's FULL site list -------
        tped = f"{a.tped_dir}/{a.target}.{a.chrom}.tped"
        pref = f"{tmp}/sel"
        C.run_selscan(tped, pref, "nsl")
        C.run_selscan(tped, pref, "ihh12")
        n_before = len(table)
        table = table.merge(C.read_nsl(f"{pref}.nsl.out")[["ID", "nsl"]],
                            on="ID", how="left")
        table = table.merge(C.read_ihh12(f"{pref}.ihh12.out")[["ID", "ihh12"]],
                            on="ID", how="left")
        if len(table) != n_before:
            sys.exit(f"ERROR: selscan merge changed the row count "
                     f"{n_before} -> {len(table)}; duplicate IDs")
        log(f"  selscan done ({time.time() - t0:.0f} s)")

        ihs_out = f"{tmp}/target.binom.ihs.out"
        C.run_ihsbin(f"{a.hmap_dir}/{a.target}.{a.chrom}.hap",
                     f"{a.hmap_dir}/{a.target}.{a.chrom}.map",
                     ihs_out, flavours["binom"])
        d = C.read_ihs_hapbin(ihs_out)[["ID", "ihs", "delihh"]]
        table = table.merge(d.rename(columns={"ihs": "ihs_binom",
                                              "delihh": "delihh_binom"}),
                            on="ID", how="left")
        if len(table) != n_before:
            sys.exit(f"ERROR: ihsbin merge changed the row count "
                     f"{n_before} -> {len(table)}; duplicate IDs")
        log(f"  ihsbin done ({time.time() - t0:.0f} s)")

        # ---- XP-EHH, one run per outgroup ------------------------------------
        xp_raw = {}
        for alt in alt_pops:
            o = run_xpehh(a.target, alt, a.chrom, a.hmap_dir, tmp, f"{tmp}/xp")
            xp, n_mono = C.xpehh_on_target_daf(
                C.read_xpehh(o), parsed[a.target][0], parsed[a.target][2])
            xp_raw[f"{a.target}_vs_{alt}"] = xp
            log(f"  xpehh {a.target} vs {alt}: {len(xp):,} sites")
            if a.keep_raw_xpehh:
                os.replace(o, f"{a.out_dir}/{a.target}_vs_{alt}.{a.chrom}.xpehh.out")
            else:
                os.remove(o)

        # ---- site list, applied LAST -----------------------------------------
        # Never before: thinning the grid would change the EHH integral.
        #
        # Keep a position when at least one of the COMPARED populations is
        # polymorphic there -- the target and its three outgroups. This is the
        # cosi2 analogue in both directions: it excludes variation private to
        # populations that are never compared (which cosi2 cannot produce) and
        # keeps the target's own private variation (which cosi2 does produce,
        # the target being in the sampled set).
        keep = np.zeros(len(table), dtype=bool)
        for p in pops:
            d = sub[p][2]
            keep |= (d > 0) & (d < 1)
        log(f"  site list (target+3 outgroup): {int(keep.sum()):,} of "
            f"{len(table):,} kept ({len(table) - int(keep.sum()):,} invariant "
            f"in all of {', '.join(pops)})")
        table = table[keep].reset_index(drop=True)

        # ---- persist ----------------------------------------------------------
        table.to_csv(f"{a.out_dir}/{a.target}.{a.chrom}.raw.tsv",
                     sep="\t", index=False, na_rep="NA")
        kv = {}
        for col in ("nsl", "ihh12", "ihs_binom", "delihh_binom"):
            v = table[col].to_numpy(dtype=float)
            ok = np.isfinite(v)
            kv[f"{col}__k"] = C.count_from_daf(table["daf"].to_numpy()[ok], n_hap)
            kv[f"{col}__v"] = v[ok]
        ids_keep = set(table["ID"])
        for pair, xp in xp_raw.items():
            s = xp[xp["ID"].isin(ids_keep)]
            tf = f"{a.out_dir}/{a.target}.{a.chrom}.xp_{pair}.tsv"
            s.to_csv(tf + ".tmp", sep="\t", index=False, na_rep="NA")
            os.replace(tf + ".tmp", tf)
            v = s["xpehh"].to_numpy(dtype=float)
            ok = np.isfinite(v)
            kv[f"xpehh_{pair}__k"] = C.count_from_daf(s["daf"].to_numpy()[ok], n_hap)
            kv[f"xpehh_{pair}__v"] = v[ok]
        tmpf = f"{a.out_dir}/{a.target}.{a.chrom}.kvals.npz.tmp"
        with open(tmpf, "wb") as fh:
            np.savez_compressed(fh, n_hap=n_hap, **kv)
        os.replace(tmpf, f"{a.out_dir}/{a.target}.{a.chrom}.kvals.npz")

        info = dict(target=a.target, chrom=a.chrom, n_hap=n_hap,
                    n_common=int(len(common)), n_kept=int(len(table)),
                    n_target_sites=int(len(parsed[a.target][1])),
                    pairs=sorted(xp_raw), outgroups=alt_pops,
                    secs=round(time.time() - t0, 1),
                    na={c: int(table[c].isna().sum()) for c in
                        ("mean_fst", "deldaf", "nsl", "ihh12",
                         "ihs_binom", "delihh_binom")},
                    tools=C.tool_versions(), job=os.environ.get("SLURM_JOB_ID"))
        with open(state + ".tmp", "w") as fh:
            json.dump(info, fh, indent=1, default=str)
        os.replace(state + ".tmp", state)
        log(f"{a.target} chr{a.chrom}: {len(table):,} sites in {info['secs']} s")
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
STEP 1 - Build a reusable neutral normalization reference.

Generates N neutral cosi2 simulations for one demographic model, computes every
component statistic on them, and writes mean/std tables with one bin per
exact derived allele count k (specific to n_hap).

TWO-PHASE so a walltime kill does not lose the whole run:
  phase A  each worker writes its per-sim (daf, value) arrays to pools/<stat>/
  phase B  a cheap pass pools those files and writes the bin CSVs

Phase B can be re-run alone with --pool-only, and more sims can be added later
with --resume without regenerating what already exists.

iHS and delihh use ihsbin --binom in production; --ihs-default adds the
non-binom flavour (finished archive comparison only; not used in production).
Each flavour gets its own bin table: a z-score is only valid against a null
built with the same tool flags. delihh is SIGNED. XP-EHH is binned on target
DAF and target-monomorphic sites are KEPT (k = 0..n_hap), per
cms_stats.XPEHH_KEEP_TARGET_MONOMORPHIC; MANIFEST records it.

Also runs a SPLIT-HALF diagnostic: bins built from each half of the sims
separately, so you can see whether N sims estimates the null stably.

Output
------
  bins/{model}_{Mb}mb_pop{target}/nsl_bin.csv
                                 /ihh12_bin.csv
                                 /ihs_binom_bin.csv
                                 /delihh_binom_bin.csv
                                 /xpehh_{target}_vs_{alt}_bin.csv
                                 /splithalf.tsv
                                 /MANIFEST.txt
  (ihs_default_bin.csv / delihh_default_bin.csv only with --ihs-default)
"""

import argparse
import glob
import os
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

import numpy as np
import pandas as pd

import cms_stats as C


# A retry must draw a DIFFERENT map window, so each attempt offsets the seed
# by this stride. Large enough never to collide with another combination's
# seed block in diag/combos.tsv (those are spaced by 10,000).
SEED_RETRY_STRIDE = 100_000_000

# make_one_sim return codes
SKIPPED = "skipped"


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


# ==============================================================================
# PHASE A - one neutral replicate
# ==============================================================================
def make_one_sim(task):
    """
    Simulate one neutral replicate, compute every statistic, and write its
    (daf, value) arrays to pools/<stat>/<sim_id>.npz.

    Returns (ok, sim_id, message).
    """
    (sim_id, par_file, pop_ids, target_pop, work_root, pools_root, seed,
     flavours, attempt_timeout, budget) = task

    # absolute: coalescent runs with cwd=sim_dir, so a relative --tped prefix
    # would land in sim_dir/<work_root>/... where nothing looks for it
    sim_dir = os.path.abspath(os.path.join(work_root, f"sim{sim_id:06d}"))
    os.makedirs(sim_dir, exist_ok=True)

    try:
        # ---- 1. simulate. Neutral = base .par with NO sweep line appended ---
        #
        # cosi2's runtime is heavy-tailed: under --genmapRandomRegions each
        # seed draws a different window of the recombination map, and ~2% of
        # vitti_gamma windows never finish at all (median 6 s, but some run
        # for hours while RSS climbs past 20 GB). A worker that draws one is
        # consumed permanently, which is what stalled arrays 13305426 and
        # 13365140. So: cap each attempt, and RESEED on retry.
        #
        # Reseeding, not repeating: the same seed draws the same window and
        # would hang identically. seed + k*SEED_RETRY_STRIDE is deterministic,
        # so the set stays reproducible - but the seed actually used is no
        # longer seed0+i, and is therefore recorded per replicate below.
        #
        # subprocess.run(timeout=) kills the child itself. (The shell `timeout`
        # used in diag/ returns 124 instead, which under `set -e` silently
        # kills the calling script - that bug ate two diagnostic sweeps.)
        prefix = os.path.join(sim_dir, "raw")
        env = dict(os.environ, COSI_NEWSIM="1", COSI_MAXATTEMPTS="5000")
        t_start = time.time()
        attempts, seed_used, timings = 0, None, []
        while True:
            spent = time.time() - t_start
            if spent >= budget:
                return (SKIPPED, sim_id,
                        f"budget {budget}s exhausted after {attempts} attempts "
                        f"(seeds tried: {timings})")
            this_seed = seed + attempts * SEED_RETRY_STRIDE
            cap = min(attempt_timeout, budget - spent)
            a0 = time.time()
            attempts += 1
            try:
                subprocess.run(
                    ["coalescent", "-p", par_file, "-v", "-g",
                     "--genmapRandomRegions", "--nTriesPerSim", "5000",
                     "--drop-singletons", ".25",
                     "--tped", prefix, "-n", "1", "-M", "-r", str(this_seed)],
                    check=True, capture_output=True, env=env, cwd=sim_dir,
                    timeout=cap,
                )
                seed_used = this_seed
                timings.append((this_seed, round(time.time() - a0, 1), "ok"))
                break
            except subprocess.TimeoutExpired:
                timings.append((this_seed, round(time.time() - a0, 1), "timeout"))
                for f in glob.glob(f"{prefix}_*.tped"):
                    os.remove(f)
                continue

        # ---- 2. flip to derived = 1 -----------------------------------------
        tpeds = {}
        for p in pop_ids:
            src = f"{prefix}_0_{p}.tped"
            if not os.path.exists(src):
                return False, sim_id, f"cosi2 produced no tped for pop {p}"
            dst = os.path.join(sim_dir, f"flip_0_{p}.tped")
            C.flip_tped(src, dst)
            os.remove(src)
            tpeds[p] = dst

        alt_pops = [p for p in pop_ids if p != target_pop]
        tgt = tpeds[target_pop]

        # ---- 3. hap/map, one pass per TPED ----------------------------------
        # parse_and_write_hap_map replaces parse_tped + write_hap_map, which
        # read the target twice. It also makes every population's positions
        # available, so 01 can assert the shared site list that 02 gets from
        # sfs_table: ONE map (the target's) is handed to every XP-EHH pair.
        mapf = os.path.join(sim_dir, "target.map")
        parsed, haps = {}, {}
        for p in pop_ids:
            h = os.path.join(sim_dir, f"pop{p}.hap")
            parsed[p] = C.parse_and_write_hap_map(
                tpeds[p], h, mapf if p == target_pop else None)
            haps[p] = h
        tgt_ids, tgt_pos, tgt_daf, n_hap = parsed[target_pop]
        for p in alt_pops:
            if not np.array_equal(parsed[p][1], tgt_pos):
                raise C.InputError(
                    f"sim {sim_id}: pop {p} position list differs from target "
                    f"{target_pop}; the target's map is used for every pair")
        C.check_unique_ids(tgt_ids, f"sim {sim_id} target TPED")
        C.check_polarization(tgt_daf, f"sim {sim_id} target pop {target_pop}")
        payload = {}

        # ---- 4. selscan: nSL, iHH12 -----------------------------------------
        pref = os.path.join(sim_dir, "sel")
        C.run_selscan(tgt, pref, "nsl")
        C.run_selscan(tgt, pref, "ihh12")
        nsl = C.read_nsl(f"{pref}.nsl.out")
        ihh12 = C.read_ihh12(f"{pref}.ihh12.out")
        # The pool stores the TOOL's frequency, so the bin tables are keyed on
        # k = rint(tool_freq * n_hap) while 02 normalizes on the TPED's DAF.
        # Asserting they agree here is what makes the two routes one route.
        for lbl, t in (("selscan --nsl", nsl), ("selscan --ihh12", ihh12)):
            C.check_unique_ids(t["ID"], f"sim {sim_id} {lbl} output")
            C.check_tool_freq(t["ID"], t["daf"], tgt_ids, tgt_daf, n_hap,
                              f"sim {sim_id} {lbl}")
        payload["nsl"] = (nsl["daf"].to_numpy(), nsl["nsl"].to_numpy())
        payload["ihh12"] = (ihh12["daf"].to_numpy(), ihh12["ihh12"].to_numpy())

        # ---- 5. iHS + signed delihh, one run per flavour ---------------------
        for flav, flags in flavours.items():
            out = os.path.join(sim_dir, f"target.{flav}.ihs.out")
            C.run_ihsbin(haps[target_pop], mapf, out, flags)
            d = C.read_ihs_hapbin(out)
            C.check_unique_ids(d["ID"], f"sim {sim_id} ihsbin {flav} output")
            C.check_tool_freq(d["ID"], d["daf"], tgt_ids, tgt_daf, n_hap,
                              f"sim {sim_id} ihsbin {flav}")
            payload[f"ihs_{flav}"] = (d["daf"].to_numpy(), d["ihs"].to_numpy())
            payload[f"delihh_{flav}"] = (d["daf"].to_numpy(), d["delihh"].to_numpy())

        # ---- 6. XP-EHH, one per pair ----------------------------------------
        for p in alt_pops:
            o = os.path.join(sim_dir, f"{target_pop}_vs_{p}.xpehh.out")
            C.run_xpehhbin(haps[target_pop], haps[p], mapf, o)
            # binned on TARGET DAF; target-monomorphic sites per C.XPEHH_KEEP_TARGET_MONOMORPHIC
            d, _ = C.xpehh_on_target_daf(C.read_xpehh(o), tgt_ids, tgt_daf)
            payload[f"xpehh_{target_pop}_vs_{p}"] = (d["daf"].to_numpy(),
                                                     d["xpehh"].to_numpy())

        # ---- 7. persist to the pool (phase B reads these) -------------------
        for stat, (daf, val) in payload.items():
            d = os.path.join(pools_root, stat)
            os.makedirs(d, exist_ok=True)
            tmp = os.path.join(d, f"{sim_id:06d}.npz.tmp")
            # file handle, not path: given a path, numpy appends ".npz" and
            # the tmp file becomes *.npz.tmp.npz, which the pool glob picks up
            with open(tmp, "wb") as fh:
                np.savez_compressed(fh, daf=daf, val=val, n_hap=n_hap)
            os.replace(tmp, os.path.join(d, f"{sim_id:06d}.npz"))

        # The seed is no longer derivable as seed0+i once a retry happens, so
        # record what was actually used. Without this the MANIFEST's
        # reproducibility claim is false for every reseeded replicate.
        meta_dir = os.path.join(pools_root, "meta")
        os.makedirs(meta_dir, exist_ok=True)
        tmp = os.path.join(meta_dir, f"{sim_id:06d}.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(dict(sim_id=sim_id, seed_base=seed, seed_used=seed_used,
                           attempts=attempts, attempt_log=timings,
                           total_s=round(time.time() - t_start, 1)), fh)
        os.replace(tmp, os.path.join(meta_dir, f"{sim_id:06d}.json"))

        return True, sim_id, f"ok seed={seed_used} attempts={attempts}"

    except C.InputError as e:
        # A guard failure fails THIS replicate; it is logged and excluded from
        # the pool rather than aborting a 7,000-replicate run.
        return False, sim_id, f"GUARD: {e}"
    except subprocess.CalledProcessError as e:
        # cosi2 prints per-attempt progress first; the actual error is at the end
        err = (e.stderr or b"").decode()[-600:]
        return False, sim_id, f"{e.cmd[0]} failed: {err}"
    except Exception as e:
        return False, sim_id, f"{type(e).__name__}: {e}"
    finally:
        shutil.rmtree(sim_dir, ignore_errors=True)


def resolve_par(par_file, region_mb, out_dir):
    """Check the .par against --region-mb and make recomb_file absolute.

    cosi2 resolves recomb_file against the cwd, and each sim runs in its own
    sim_dir, so a relative path fails on every attempt. Writes the resolved
    copy into out_dir, which also records exactly what was simulated.
    """
    par_dir = os.path.dirname(os.path.abspath(par_file))
    lines, length = [], None
    with open(par_file) as f:
        for line in f:
            tok = line.split()
            if tok and tok[0] == "length":
                length = int(tok[1])
            if tok and tok[0] == "recomb_file":
                recom = tok[1] if os.path.isabs(tok[1]) else os.path.join(par_dir, tok[1])
                if not os.path.exists(recom):
                    sys.exit(f"recomb_file not found: {recom}")
                line = f"recomb_file {recom}\n"
            lines.append(line)
    if length != region_mb * 1_000_000:
        sys.exit(f"{par_file} has length {length}, but --region-mb is {region_mb}")
    out = os.path.join(out_dir, "resolved.par")
    with open(out, "w") as f:
        f.writelines(lines)
    return os.path.abspath(out)


def already_done(pools_root, stat_hint, sim_id):
    return os.path.exists(os.path.join(pools_root, stat_hint,
                                       f"{sim_id:06d}.npz"))


# ==============================================================================
# PHASE B - pool and bin
# ==============================================================================
def load_pool(pools_root, stat, sim_subset=None):
    """Concatenate every per-sim npz for one statistic."""
    files = sorted(glob.glob(os.path.join(pools_root, stat, "*.npz")))
    if sim_subset is not None:
        keep = {f"{s:06d}.npz" for s in sim_subset}
        files = [f for f in files if os.path.basename(f) in keep]
    dafs, vals, nhaps = [], [], set()
    for f in files:
        try:
            with np.load(f) as z:
                dafs.append(z["daf"])
                vals.append(z["val"])
                if "n_hap" in z:
                    nhaps.add(int(z["n_hap"]))
        except Exception as e:
            log(f"  WARNING: unreadable pool file {f}: {e}")
    if not dafs:
        return None, None, 0, None
    if len(nhaps) > 1:
        raise RuntimeError(f"{pools_root}/{stat}: pool mixes n_hap {sorted(nhaps)}; "
                           f"counts are not comparable across sample sizes")
    return (np.concatenate(dafs), np.concatenate(vals), len(files),
            nhaps.pop() if nhaps else None)


def pool_and_bin(pools_root, out_dir, n_hap_arg=None):
    """Phase B. Returns a summary list and the set of statistics found."""
    stats = sorted(d for d in os.listdir(pools_root)
                   if os.path.isdir(os.path.join(pools_root, d)))
    if not stats:
        raise RuntimeError(f"no pooled data in {pools_root}")

    summary = []
    n_haps = set()
    for stat in stats:
        daf, val, nsim, nh = load_pool(pools_root, stat)
        if daf is None:
            log(f"  WARNING: no data for {stat}")
            continue
        nh = nh if nh is not None else n_hap_arg
        if nh is None:
            raise RuntimeError(
                f"{stat}: pool files carry no n_hap (written before exact-count "
                f"binning) and --n-hap was not given. Counts cannot be derived "
                f"from DAF without it.")
        if n_hap_arg is not None and nh != n_hap_arg:
            raise RuntimeError(f"{stat}: pool says n_hap={nh}, --n-hap says {n_hap_arg}")
        n_haps.add(nh)
        out_csv = os.path.join(out_dir, f"{stat}_bin.csv")
        tbl = C.build_bins(daf, val, nh, out_csv)
        if tbl is None:
            log(f"  WARNING: {stat} produced no usable bins")
            continue
        usable = int(tbl["mean"].notna().sum())
        thin = int(tbl.loc[tbl["mean"].notna(), "count"].min()) if usable else 0
        log(f"  {stat:<26} {nsim:>5} sims | {len(daf):>11,} sites | "
            f"{usable}/{len(tbl)} usable counts | thinnest {thin:,}")
        summary.append((stat, nsim, len(daf), usable, thin))
    if len(n_haps) > 1:
        raise RuntimeError(f"statistics disagree on n_hap: {sorted(n_haps)}")
    return summary, stats, (n_haps.pop() if n_haps else None)


def split_half(pools_root, out_dir, stats, n_hap):
    """
    Build bins from each half of the sims and compare, to measure whether the
    sim count estimates the null stably. Large half-to-half differences mean
    the null is still noisy and more sims are needed.
    """
    rows = []
    for stat in stats:
        files = sorted(glob.glob(os.path.join(pools_root, stat, "*.npz")))
        ids = sorted(int(os.path.basename(f).split(".")[0]) for f in files)
        if len(ids) < 4:
            continue
        a, b = ids[0::2], ids[1::2]          # interleaved, not first/second half

        tabs = []
        for subset in (a, b):
            daf, val, _, _ = load_pool(pools_root, stat, subset)
            if daf is None:
                tabs.append(None)
                continue
            tabs.append(C.build_bins(daf, val, n_hap))
        if any(t is None for t in tabs):
            continue

        m = tabs[0].merge(tabs[1], on="k", suffixes=("_a", "_b"))
        m = m[m["mean_a"].notna() & m["mean_b"].notna()]
        if m.empty:
            continue

        # Difference in bin mean, scaled by the bin sigma: how far a z-score
        # would move if you had drawn the other half as your null.
        shift = ((m["mean_a"] - m["mean_b"]).abs()
                 / m[["std_a", "std_b"]].mean(axis=1))
        sd_ratio = (m["std_a"] / m["std_b"])

        rows.append(dict(
            statistic=stat,
            bins_compared=len(m),
            max_mean_shift_in_sd=round(float(shift.max()), 4),
            median_mean_shift_in_sd=round(float(shift.median()), 4),
            min_sd_ratio=round(float(sd_ratio.min()), 4),
            max_sd_ratio=round(float(sd_ratio.max()), 4),
        ))

    if not rows:
        return None
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_dir, "splithalf.tsv"), sep="\t", index=False)

    log("")
    log("split-half stability (two independent halves of the neutral set):")
    log(df.to_string(index=False))
    worst = df["max_mean_shift_in_sd"].max()
    log("")
    log(f"  worst bin-mean shift between halves: {worst:.3f} sigma")
    if worst > 0.25:
        log("  -> The null is still noisy at this sim count. Z-scores would move")
        log("     materially if you had drawn a different neutral set. Consider")
        log("     more sims before committing to these bins.")
    else:
        log("  -> The null looks stable at this sim count.")
    return df


# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Generate neutral sims and build normalization bins.")
    ap.add_argument("--model", required=True, help="e.g. gravel")
    ap.add_argument("--par-file", required=True,
                    help="cosi2 .par for this model, with NO sweep line")
    ap.add_argument("--target-pop", required=True)
    ap.add_argument("--pop-ids", required=True, help="comma separated, e.g. 1,2,3")
    ap.add_argument("--region-mb", type=int, required=True,
                    help="region length in Mb; must match the .par file")
    ap.add_argument("--n-sims", type=int, default=800)
    ap.add_argument("--bins-root", default="bins")
    ap.add_argument("--pools-root", default=None,
                    help="where per-sim arrays live (default: pools/<bins dir name>)")
    ap.add_argument("--work-root", default="scratch_neutral")
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--seed0", type=int, default=0,
                    help="base RNG seed; sim i uses seed0+i (0 = cosi2 default)")
    ap.add_argument("--resume", action="store_true",
                    help="skip sims already present in the pool")
    ap.add_argument("--n-hap", type=int, default=None,
                    help="haplotypes per population; only needed to rebuild a "
                         "pool written before n_hap was stored in the npz")
    ap.add_argument("--ihs-default", action="store_true",
                    help="also run the non-binom iHS/delihh flavour (archive "
                         "comparison only; not used in production)")
    ap.add_argument("--sim-timeout", type=float, default=120.0,
                    help="seconds per cosi2 ATTEMPT before it is killed and "
                         "retried with a fresh seed (default 120; the median "
                         "replicate is ~6 s and the slowest finisher "
                         "measured was 11 s)")
    ap.add_argument("--sim-budget", type=float, default=600.0,
                    help="total seconds per replicate across all attempts; "
                         "once exhausted the replicate is SKIPPED and logged "
                         "to pools/<combo>/skipped.tsv (default 600)")
    ap.add_argument("--pool-only", action="store_true",
                    help="skip simulation; just rebuild bins from the pool")
    args = ap.parse_args()
    flavours = C.ihs_flavours(args.ihs_default)

    pop_ids = [p.strip() for p in args.pop_ids.split(",")]
    if args.target_pop not in pop_ids:
        sys.exit(f"target pop {args.target_pop} not in pop_ids {pop_ids}")

    out_dir = C.bins_dir_for(args.bins_root, args.model,
                             args.region_mb, args.target_pop)
    os.makedirs(out_dir, exist_ok=True)

    pools_root = args.pools_root or os.path.join("pools", os.path.basename(out_dir))
    os.makedirs(pools_root, exist_ok=True)

    log(f"model {args.model} | {args.region_mb} Mb | target pop {args.target_pop} "
        f"| pops {pop_ids}")
    log(f"bins  -> {out_dir}")
    log(f"pools -> {pools_root}")

    # ---------------- PHASE A ------------------------------------------------
    if not args.pool_only:
        if not os.path.exists(args.par_file):
            sys.exit(f"par file not found: {args.par_file}")
        try:
            C.check_tools(["coalescent", "selscan", "ihsbin", "xpehhbin"])
        except RuntimeError as e:
            sys.exit(f"ERROR: {e}")

        os.makedirs(args.work_root, exist_ok=True)
        par_file = resolve_par(args.par_file, args.region_mb, out_dir)

        tasks = []
        skipped = 0
        for i in range(args.n_sims):
            if args.resume and already_done(pools_root, "nsl", i):
                skipped += 1
                continue
            tasks.append((i, par_file, pop_ids,
                          args.target_pop, args.work_root, pools_root,
                          args.seed0 + i if args.seed0 else 0, flavours,
                          args.sim_timeout, args.sim_budget))

        if skipped:
            log(f"resume: {skipped} sims already pooled, {len(tasks)} to run")
        log(f"simulating {len(tasks)} sims with {args.workers} workers "
            f"(iHS flavours: {list(flavours)})")

        ok, failed, skipped = 0, [], []
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(make_one_sim, t) for t in tasks]
            for fut in as_completed(futs):
                good, sid, msg = fut.result()
                done += 1
                if good is True:
                    ok += 1
                elif good == SKIPPED:
                    skipped.append((sid, msg))
                else:
                    failed.append((sid, msg))
                # progress must count EVERY outcome, not just successes:
                # a run where nothing succeeds otherwise logs nothing at all
                # for its whole wall time, which is how arrays 13305426 and
                # 13365140 looked healthy while producing zero replicates
                if done % 50 == 0:
                    log(f"  {done}/{len(tasks)} attempted "
                        f"({ok} ok, {len(skipped)} skipped, {len(failed)} failed)")

        log(f"simulation complete: {ok} ok, {len(skipped)} skipped, "
            f"{len(failed)} failed")
        for sid, msg in failed[:10]:
            log(f"  FAILED sim {sid}: {msg}")

        # The skip log is the realized blow-up rate, which is the input to
        # the question of whether excluding those replicates biases the null.
        if skipped:
            skip_path = os.path.join(pools_root, "skipped.tsv")
            new = not os.path.exists(skip_path)
            with open(skip_path, "a") as f:
                if new:
                    f.write("sim_id\tdetail\n")
                for sid, msg in skipped:
                    f.write(f"{sid}\t{msg}\n")
            log(f"  {len(skipped)} skipped replicates logged to {skip_path}")
            for sid, msg in skipped[:5]:
                log(f"  SKIPPED sim {sid}: {msg}")
        shutil.rmtree(args.work_root, ignore_errors=True)

    # ---------------- PHASE B ------------------------------------------------
    log("pooling and binning ...")
    summary, stats, n_hap = pool_and_bin(pools_root, out_dir, args.n_hap)

    log("")
    sh = split_half(pools_root, out_dir, stats, n_hap)

    # ---------------- MANIFEST ------------------------------------------------
    n_pooled = max((s[1] for s in summary), default=0)
    with open(os.path.join(out_dir, "MANIFEST.txt"), "w") as f:
        f.write(f"generated      : {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        f.write(f"model          : {args.model}\n")
        f.write(f"region_mb      : {args.region_mb}\n")
        f.write(f"target_pop     : {args.target_pop}\n")
        f.write(f"pop_ids        : {','.join(pop_ids)}\n")
        f.write(f"par_file       : {os.path.abspath(args.par_file)}\n")
        f.write(f"n_sims_asked   : {args.n_sims}\n")
        f.write(f"n_sims_pooled  : {n_pooled}\n")
        f.write(f"seed0          : {args.seed0} "
                f"({'sim i uses seed0+i' if args.seed0 else 'unseeded, cosi2 random'})\n")
        f.write(f"sim_timeout    : {args.sim_timeout}s per cosi2 attempt, "
                f"then killed and retried at seed + "
                f"attempt*{SEED_RETRY_STRIDE}\n")
        f.write(f"sim_budget     : {args.sim_budget}s per replicate, then "
                f"SKIPPED\n")
        # seed0+i holds only for replicates that never needed a retry; the
        # seed actually used is in pools/<combo>/meta/<id>.json
        n_meta = len(glob.glob(os.path.join(pools_root, "meta", "*.json")))
        n_retried = 0
        for mf in glob.glob(os.path.join(pools_root, "meta", "*.json")):
            try:
                with open(mf) as mh:
                    if json.load(mh).get("attempts", 1) > 1:
                        n_retried += 1
            except Exception:                                # noqa: BLE001
                pass
        f.write(f"seeds recorded : {n_meta} replicates have meta/*.json; "
                f"{n_retried} needed a reseed (seed0+i is NOT the seed used "
                f"for those)\n")
        skip_path = os.path.join(pools_root, "skipped.tsv")
        n_skip = (sum(1 for _ in open(skip_path)) - 1) if os.path.exists(skip_path) else 0
        f.write(f"skipped        : {n_skip} replicates exhausted the budget "
                f"(see skipped.tsv)\n")
        f.write(f"n_hap          : {n_hap}\n")
        f.write(f"bin key        : exact derived count k, 0..{n_hap}\n")
        f.write(f"min_per_bin    : {C.MIN_SITES_PER_BIN}\n")
        f.write("neutral        : base .par, NO sweep line appended\n")
        f.write("cosi2 flags    : --genmapRandomRegions --nTriesPerSim 5000 "
                "--drop-singletons .25 -n 1 -M\n")
        f.write("coding         : flipped to derived = 1\n")
        f.write("delihh         : SIGNED (iHH_1 - iHH_0)\n")
        vers = C.tool_versions()
        f.write(f"selscan        : {vers['selscan']}\n")
        f.write(f"hapbin         : {vers['hapbin']}\n")
        f.write(f"nsl/ihh12 cmd  : selscan --nsl|--ihh12 --threads 1 {' '.join(C.SELSCAN_FLAGS)}\n")
        for flav, flags in flavours.items():
            f.write(f"ihs flavour    : {flav} -> ihsbin {' '.join([*C.IHSBIN_FLAGS, *flags])}\n")
        f.write(f"xpehh cmd      : xpehhbin {' '.join(C.XPEHHBIN_FLAGS)}\n")
        if C.XPEHH_KEEP_TARGET_MONOMORPHIC:
            f.write("xpehh bin freq : TARGET pop DAF (joined on ID); k = 0..n_hap\n")
            f.write("XPEHH_KEEP_TARGET_MONOMORPHIC : True\n")
        else:
            f.write("xpehh bin freq : TARGET pop DAF (joined on ID); sites monomorphic "
                    "in the target dropped\n")
        f.write("binning rule   : one bin per exact derived count k (no frequency\n"
                "                 intervals). Equal-width bins left up to ~0.95 z of\n"
                "                 residual frequency dependence INSIDE a bin.\n")
        f.write("\nstatistic, n_sims, n_sites, usable_bins, thinnest_usable_bin\n")
        for row in summary:
            f.write(", ".join(str(x) for x in row) + "\n")
        if sh is not None:
            f.write("\nsplit-half stability: see splithalf.tsv\n")
            f.write(f"worst bin-mean shift between halves (sigma): "
                    f"{sh['max_mean_shift_in_sd'].max():.4f}\n")

    log("")
    log(f"done -> {out_dir}")
    log(f"pool retained at {pools_root} (delete when the bins are accepted, or "
        f"keep to extend the set later with --resume)")


if __name__ == "__main__":
    main()

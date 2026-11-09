#!/usr/bin/env python3
"""
STEP 2 - Recompute every component statistic for one TPED batch (archived,
replacement or refill), normalized against the shared neutral bins from step 1.

Fixes carried in, relative to the original pipeline:
  * deldaf    = D_target - mean(D_outgroups)              (was max - min)
  * max_xpehh joined on variant ID per site, not row index (was index-aligned)
  * FST       = allele-count variance-component estimator (was Nei G_ST);
                fixed-same pairs are NaN and excluded from mean_fst (CMS)
  * nSL output header auto-detected                       (was skiprows=1 always)
  * one shared neutral null instead of a per-batch draw

iHS and delihh: production is *_binom only (ihsbin --binom, which removes
the 1/k EHH floor and widens coverage from DAF 0.105-0.895 to 0.05-0.95).
--ihs-default adds the *_default flavour; that was for the finished
null-to-null archive comparison and is not used in production.
delihh is SIGNED. XP-EHH is binned on target DAF with target-monomorphic
sites KEPT; every pair is required.

Hard checks (any failure aborts the batch; nothing is written): TPED
identity (zip/metafile names, member paths, per-sim pop set, no identical
pops), target label + causal DAF gate [0.05, 0.95], one metafile row per
sim, no failed or unwritten sim, every XP-EHH pair present.
"""

import argparse
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

import numpy as np
import pandas as pd

import cms_stats as C

GCS_BASE = "gs://deepsweep/raw_sim_tped_2026"
TPED_RE = re.compile(r"sel\.hap\.(\d+)_0_(\d+)\.tped")
TPED_MEMBER_RE = re.compile(r"sel/sel\.hap\.(\d+)_0_(\d+)\.tped")   # whole path

# Sims with NO causal row anywhere (verified over all 1,028 archived batches:
# pos == region/2 holds for 1,016,901/1,016,901 causal rows, and exactly these
# two sims have none). Dropped by decision 2026-10-15.
DROP_SIMS = {
    ("vitti_beta", "13050001012"): {"715"},
    ("vitti_gamma", "4501504016"): {"517"},
}


def out_cols(flavours):
    cols = ["sim_batch_no", "sim_id", "pos",
            "mean_fst", "deldaf", "daf", "maf",
            "norm_nsl", "norm_ihh12"]
    for flav in flavours:
        cols += [f"norm_ihs_{flav}", f"norm_delihh_{flav}"]
    return cols + ["max_xpehh", "deri_gen", "sel_gen", "s", "selpos"]


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True)


# ==============================================================================
# PER-SIM WORKER
# ==============================================================================
LABEL_FAIL = "LABEL"      # prefix of a failure that aborts the whole batch
IDENTITY_FAIL = "IDENTITY"   # TPED is not the (sim, pop) that was asked for
GUARD_FAIL = "GUARD"      # a cms_stats.InputError: a statistic would be wrong


def component_cols(flavours):
    """The nine components -- everything in the output that is not a key."""
    keys = {"sim_batch_no", "sim_id", "pos", "deri_gen", "sel_gen", "s",
            "selpos"}
    return [c for c in out_cols(flavours) if c not in keys]
CAUSAL_DAF_MIN, CAUSAL_DAF_MAX = 0.05, 0.95   # the corpus window, inclusive


def check_target_label(parsed, target_pop, selpos, expected_daf):
    """Is the sweep really in the declared target population?

    The batch id declares the target; ten archived gravel batches
    (23021002001-010) declared pop 2 with the sweep in pop 1 (CLAUDE.md,
    2026-10-17). Returns None if fine, else a LABEL_FAIL message.
      1. the target has exactly one site at selpos (the causal variant)
      2. causal MAF gate: target DAF at selpos within [0.05, 0.95]
      3. the metafile's causal DAF, when it has one, equals the target's.
         The old pipeline wrote it from the pop the sweep was simulated in.
    On failure, names any other pop whose DAF at selpos fits instead.
    """
    def daf_at(p):
        hit = np.flatnonzero(np.asarray(parsed[p][1]) == selpos)
        return float(parsed[p][2][hit[0]]) if len(hit) == 1 else None, len(hit)

    d, n = daf_at(target_pop)
    others = {p: daf_at(p)[0] for p in parsed if p != target_pop}
    where = " ".join(f"pop{p}={v:.3f}" if v is not None else f"pop{p}=NA"
                     for p, v in sorted(others.items()))
    if n != 1:
        return (f"{LABEL_FAIL}: declared target pop {target_pop} has {n} sites at "
                f"selpos {selpos} (need exactly 1); other pops: {where}")
    if not (CAUSAL_DAF_MIN <= round(d, 6) <= CAUSAL_DAF_MAX):
        return (f"{LABEL_FAIL}: causal DAF {d:.3f} in declared target pop "
                f"{target_pop} is outside [{CAUSAL_DAF_MIN}, {CAUSAL_DAF_MAX}]; "
                f"other pops: {where}")
    if expected_daf is not None and not np.isclose(d, expected_daf, atol=1e-6):
        fits = [p for p, v in others.items()
                if v is not None and np.isclose(v, expected_daf, atol=1e-6)]
        return (f"{LABEL_FAIL}: metafile causal DAF {expected_daf:.6f} != "
                f"{d:.6f} in declared target pop {target_pop}; matches pop(s) "
                f"{fits or 'none'}; other pops: {where}")
    return None


def process_sim(task):
    """Compute raw statistics for one sim; write an un-normalized table."""
    (sim_id, pop_members, target_pop, work_dir, flavours, n_hap, zip_path,
     selpos, expected_daf, ref_pops) = task

    raw_dir = os.path.join(work_dir, "rawstats")
    tmp_dir = os.path.join(work_dir, "scratch", str(sim_id))
    os.makedirs(tmp_dir, exist_ok=True)

    try:
        if target_pop not in pop_members:
            return False, sim_id, f"no tped for target pop {target_pop}"

        pop_ids = sorted(pop_members)
        alt_pops = [p for p in pop_ids if p != target_pop]
        if sorted(pop_members) != sorted(ref_pops):
            return False, sim_id, f"pop set {pop_ids} != reference {sorted(ref_pops)}"

        # ---- extract THIS sim only -------------------------------------------
        # Extracting the whole batch up front costs 18.8 GB (gravel 3 Mb) to
        # 69.1 GB (medina-munoz 5 Mb) resident for the entire run, and that
        # decompression is serial. Per-sim extraction bounds peak disk to
        # workers x one sim's TPEDs and moves the decompression into the
        # parallel workers.
        pop_files, digests = {}, {}
        with zipfile.ZipFile(zip_path) as zf:
            for p, member in pop_members.items():
                # identity: the member must be exactly (this sim, this pop)
                m = TPED_MEMBER_RE.fullmatch(member)
                if not m or (m.group(1), m.group(2)) != (str(sim_id), str(p)):
                    return False, sim_id, (f"{IDENTITY_FAIL}: asked for sim {sim_id} "
                                           f"pop {p}, member is {member}")
                dst = os.path.join(tmp_dir, f"pop{p}.tped")
                h = hashlib.md5()
                with zf.open(member) as src, open(dst, "wb") as out:
                    for chunk in iter(lambda: src.read(1 << 20), b""):
                        h.update(chunk); out.write(chunk)
                pop_files[p], digests[p] = dst, h.hexdigest()
        # Pops of one sim share positions but never genotypes: byte-identical
        # TPEDs mean one population's file was copied or mislabelled.
        seen = {}
        for p, dg in sorted(digests.items()):
            if dg in seen:
                return False, sim_id, (f"{IDENTITY_FAIL}: sim {sim_id} pops {seen[dg]} "
                                       f"and {p} have byte-identical TPEDs")
            seen[dg] = p

        # ---- one pass per TPED: parse AND write hap/map -----------------------
        mapf = os.path.join(tmp_dir, "target.map")
        parsed, haps = {}, {}
        for p in pop_ids:
            h = os.path.join(tmp_dir, f"pop{p}.hap")
            parsed[p] = C.parse_and_write_hap_map(
                pop_files[p], h, mapf if p == target_pop else None)
            haps[p] = h

        if len(parsed[target_pop][1]) == 0:
            return False, sim_id, "target tped empty"
        if parsed[target_pop][3] != n_hap:
            return False, sim_id, (f"target pop has {parsed[target_pop][3]} haplotypes "
                                   f"but the reference set was built at {n_hap}")

        # ---- review guards, 2026-10-28 --------------------------------------
        # (3) IDs are the join key for every statistic downstream.
        for p in pop_ids:
            C.check_unique_ids(parsed[p][0], f"sim {sim_id} pop {p} TPED")
        # (2) nothing else asserts the archive arrives flipped to derived = 1.
        C.check_polarization(parsed[target_pop][2],
                             f"sim {sim_id} target pop {target_pop}")
        # (4) hapbin is handed pop p's .hap with the TARGET's .map; they must
        # at least describe the same number of sites. sfs_table's position
        # assertion validates `parsed`, not these files.
        n_map = sum(1 for _ in open(mapf))
        for p in pop_ids:
            n_h = sum(1 for _ in open(haps[p]))
            if n_h != n_map or n_h != len(parsed[p][1]):
                return False, sim_id, (
                    f"{IDENTITY_FAIL}: pop {p} hap has {n_h} rows, target map "
                    f"{n_map}, parsed {len(parsed[p][1])}")

        bad = check_target_label(parsed, target_pop, selpos, expected_daf)
        if bad:
            return False, sim_id, bad

        # SFS, keyed by variant ID so the haplotype stats can join on it.
        # sfs_table asserts that every pop shares the target's position list.
        table = C.sfs_table(parsed, target_pop, alt_pops)

        # ---- selscan: nSL, iHH12 (the last consumers of the TPED) ------------
        pref = os.path.join(tmp_dir, "sel")
        C.run_selscan(pop_files[target_pop], pref, "nsl")
        C.run_selscan(pop_files[target_pop], pref, "ihh12")
        # (1) the bins are keyed on k from the TOOL's printed frequency, the
        # batch on k from the TPED's DAF. check_tool_freq asserts the two
        # routes land on the same count at every site.
        tgt_ids, tgt_daf = parsed[target_pop][0], parsed[target_pop][2]
        for stat, path, reader in (("nsl", f"{pref}.nsl.out", C.read_nsl),
                                   ("ihh12", f"{pref}.ihh12.out", C.read_ihh12)):
            t = reader(path)
            C.check_unique_ids(t["ID"], f"sim {sim_id} selscan --{stat} output")
            C.check_tool_freq(t["ID"], t["daf"], tgt_ids, tgt_daf, n_hap,
                              f"sim {sim_id} selscan --{stat}")
            table = C.safe_merge(table, t[["ID", stat]], "ID", "left",
                                 f"sim {sim_id} {stat}")

        # TPEDs are dead once selscan is done; hapbin reads only .hap/.map.
        # Dropping them here roughly halves this worker's disk footprint for
        # the whole hapbin stage, which is the longer half at 4 populations.
        for f in pop_files.values():
            try:
                os.remove(f)
            except OSError:
                pass

        # ---- iHS + signed delihh, BOTH flavours ------------------------------
        for flav, flags in flavours.items():
            out = os.path.join(tmp_dir, f"target.{flav}.ihs.out")
            C.run_ihsbin(haps[target_pop], mapf, out, flags)
            dd = C.read_ihs_hapbin(out)
            C.check_unique_ids(dd["ID"], f"sim {sim_id} ihsbin {flav} output")
            # hapbin reports the frequency of the .hap it was given, so this
            # is also the check that haps[target] really is the target pop.
            C.check_tool_freq(dd["ID"], dd["daf"], tgt_ids, tgt_daf, n_hap,
                              f"sim {sim_id} ihsbin {flav}")
            d = dd[["ID", "ihs", "delihh"]].rename(
                columns={"ihs": f"ihs_{flav}", "delihh": f"delihh_{flav}"})
            table = C.safe_merge(table, d, "ID", "left",
                                 f"sim {sim_id} ihs_{flav}")

        # ---- XP-EHH ----------------------------------------------------------
        for p in alt_pops:
            o = os.path.join(raw_dir, f"{sim_id}_{target_pop}_vs_{p}.xpehh.out")
            C.run_xpehhbin(haps[target_pop], haps[p], mapf, o)

        table.to_csv(os.path.join(raw_dir, f"{sim_id}.raw.tsv"),
                     sep="\t", index=False)
        return True, sim_id, "ok"

    except C.InputError as e:
        return False, sim_id, f"{GUARD_FAIL}: {e}"
    except subprocess.CalledProcessError as e:
        return False, sim_id, f"{e.cmd[0]} failed: {(e.stderr or b'').decode()[:200]}"
    except Exception as e:
        return False, sim_id, f"{type(e).__name__}: {e}"
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def finalize_sim(task):
    """Wrapper: a guard failure must fail THIS SIM, which makes written !=
    len(tasks) and aborts the batch, rather than raising out of the pool."""
    try:
        return _finalize_sim(task)
    except C.InputError as e:
        return 0, [f"  sim {task[0]}: FAILED {GUARD_FAIL}: {e}"]


def _finalize_sim(task):
    """Normalize one sim's raw table against the bins and write its TSV.
    Returns (1 if written else 0, [log messages])."""
    (sid, work_dir, flavours, bins, n_hap, pairs, batch_id, selpos, model,
     meta_row) = task
    msgs = []
    raw_f = os.path.join(work_dir, "rawstats", f"{sid}.raw.tsv")
    if not os.path.exists(raw_f):
        return 0, [f"  sim {sid}: FAILED no raw table"]
    df = pd.read_csv(raw_f, sep="\t", dtype={"ID": str})

    # one-pop statistics, each against its own matching bin table
    norm_map = [("nsl", "norm_nsl", "nsl"),
                ("ihh12", "norm_ihh12", "ihh12")]
    for flav in flavours:
        norm_map += [(f"ihs_{flav}", f"norm_ihs_{flav}", f"ihs_{flav}"),
                     (f"delihh_{flav}", f"norm_delihh_{flav}", f"delihh_{flav}")]

    for src, dst, binkey in norm_map:
        if src in df.columns:
            df[dst] = C.normalize(df[src].to_numpy(dtype=float),
                                  df["daf"].to_numpy(), bins[binkey], n_hap)
        else:
            df[dst] = np.nan

    # XP-EHH: per pair, joined on variant ID (per-pair outputs legitimately
    # cover different site sets, so row order is not a safe key), then the
    # row-wise max of the normalized scores.
    xp_cols = []
    for i, pair in enumerate(pairs):
        f = os.path.join(work_dir, "rawstats", f"{sid}_{pair}.xpehh.out")
        # Every pair is required: max_xpehh over fewer pairs is a different
        # statistic. NA at individual sites stays allowed.
        if not os.path.exists(f):
            return 0, [f"  sim {sid}: FAILED missing XP-EHH output for {pair}"]
        try:
            # binned on TARGET DAF, exactly as the null; target-monomorphic
            # sites kept/dropped per C.XPEHH_KEEP_TARGET_MONOMORPHIC (k = 0 / n_hap)
            d, _ = C.xpehh_on_target_daf(C.read_xpehh(f), df["ID"], df["daf"])
        except Exception as e:
            return 0, [f"  sim {sid}: FAILED XP-EHH {pair} unreadable ({e})"]
        col = f"_xp{i}"
        d[col] = C.normalize(d["xpehh"].to_numpy(), d["daf"].to_numpy(),
                             bins[f"xpehh_{pair}"], n_hap)
        df = C.safe_merge(df, d[["ID", col]], "ID", "left",
                          f"sim {sid} xpehh {pair}")
        xp_cols.append(col)

    df["max_xpehh"] = df[xp_cols].max(axis=1) if xp_cols else np.nan
    if xp_cols:
        df = df.drop(columns=xp_cols)

    # sweep parameters from the metafile, on the selected site only
    df["sim_batch_no"] = batch_id
    df["sim_id"] = sid
    df["selpos"] = (df["pos"] == selpos).astype(int)
    for c in ("deri_gen", "sel_gen", "s"):
        df[c] = np.nan
    if meta_row is not None:
        mask = df["pos"] == selpos
        for c in ("deri_gen", "sel_gen", "s"):
            if c in meta_row:
                df.loc[mask, c] = meta_row[c]

    # (5) reindex fabricates an all-NaN column for any name not in the frame,
    # so a renamed intermediate ships a silently empty component instead of
    # failing. Require every expected name to exist, and no component to be
    # entirely NaN. Measured over 5,574 real sims from all six models: 0
    # missing columns, 0 all-NaN components, so neither can fire on good data.
    want = out_cols(flavours)
    absent = [c for c in want if c not in df.columns]
    if absent:
        return 0, [f"  sim {sid}: FAILED {GUARD_FAIL}: expected column(s) "
                   f"{absent} were never produced; reindex would have "
                   f"written them as all-NaN"]
    empty = [c for c in component_cols(flavours) if df[c].isna().all()]
    if empty:
        return 0, [f"  sim {sid}: FAILED {GUARD_FAIL}: component(s) {empty} "
                   f"are entirely NaN over {len(df)} sites"]

    out = df.reindex(columns=want).round(6)
    out.to_csv(os.path.join(work_dir, "tsv",
                            f"{model}_batch{batch_id}_cms_stats_{sid}.tsv"),
               sep="\t", index=False, na_rep="NA")
    return 1, msgs


# ==============================================================================
def load_bins(bins_dir, pairs, flavours):
    """Load every bin table, failing loudly if one is missing."""
    need = ["nsl", "ihh12"]
    for flav in flavours:
        need += [f"ihs_{flav}", f"delihh_{flav}"]
    need += [f"xpehh_{p}" for p in pairs]

    bins = {}
    for key in need:
        f = os.path.join(bins_dir, f"{key}_bin.csv")
        if not os.path.exists(f):
            raise FileNotFoundError(
                f"missing bin table: {f}\n"
                f"Run 01_make_neutral_bins.py for this model/region/pop first.")
        bins[key] = pd.read_csv(f)
    return bins


def main():
    """Wrapper: guarantees scratch removal on every exit path."""
    work_dir_holder = {}
    try:
        return _main(work_dir_holder)
    finally:
        wd = work_dir_holder.get("path")
        if wd and os.path.isdir(wd):
            shutil.rmtree(wd, ignore_errors=True)
            log(f"cleaned scratch {wd}")


def _main(work_dir_holder=None):
    ap = argparse.ArgumentParser(
        description="Recompute all component stats for one archived batch.")
    ap.add_argument("--model", default="gravel")
    ap.add_argument("--batch-id", required=True)
    ap.add_argument("--bins-root", default="bins")
    ap.add_argument("--work-root", default=".")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--keep-scratch", action="store_true")
    ap.add_argument("--ihs-default", action="store_true",
                    help="also compute the non-binom iHS/delihh flavour "
                         "(archive comparison only)")
    ap.add_argument("--drop-list", default=None,
                    help="TSV with model, sim_batch_no, sim_id (e.g. "
                         "repl_sims/dropped_archive_sims.tsv): sims of THIS "
                         "(model, batch) listed there are skipped and logged")
    ap.add_argument("--input-dir", default=None,
                    help="read the tped zip + par_inputs tsv from here (made by "
                         "fetch_batch.sh on the login node) instead of GCS; "
                         "required on compute nodes, which have no internet")
    args = ap.parse_args()
    flavours = C.ihs_flavours(args.ihs_default)

    try:
        C.check_tools(["selscan", "ihsbin", "xpehhbin"])
    except RuntimeError as e:
        sys.exit(f"ERROR: {e}")

    batch_id = args.batch_id
    target_pop, region_mb, selpos = C.decode_batch_id(batch_id)
    out_dir = args.out_dir or f"{args.model}_recomputed"
    bins_dir = C.bins_dir_for(args.bins_root, args.model, region_mb, target_pop)

    log(f"batch {batch_id} | model {args.model} | target pop {target_pop} "
        f"| {region_mb} Mb | selpos {selpos:,}")
    log(f"bins: {bins_dir}")
    if not os.path.isdir(bins_dir):
        sys.exit(f"ERROR: no neutral bins at {bins_dir} - run step 1 first")
    # The MANIFEST is append-only provenance: its original section can state
    # settings that a later dated block supersedes (e.g. "target-monomorphic
    # dropped" before the 2026-10-16 XP-EHH rebuild). Echoing it verbatim put
    # superseded lines into every batch log, so log its identity (it is also
    # copied into the output zip) and the settings this run actually uses.
    _manf = os.path.join(bins_dir, "MANIFEST.txt")
    if os.path.exists(_manf):
        import hashlib as _h
        log(f"bin MANIFEST: {_manf} md5 {_h.md5(open(_manf, 'rb').read()).hexdigest()} "
            f"(copied into the output zip)")
    log("--- effective settings (this run) ---")
    log(f"  selscan {' '.join(C.SELSCAN_FLAGS)} | ihsbin {' '.join(C.IHSBIN_FLAGS)} "
        f"+ flavours {list(flavours)} | xpehhbin {' '.join(C.XPEHHBIN_FLAGS)}")
    log(f"  XP-EHH: binned on TARGET DAF, target-monomorphic sites "
        f"{'KEPT' if C.XPEHH_KEEP_TARGET_MONOMORPHIC else 'DROPPED'}; every pair required")
    log(f"  norm: one bin per exact derived count k; MIN_SITES_PER_BIN "
        f"{C.MIN_SITES_PER_BIN}")
    log("  FST: Weir-Hill, fixed-same pairs NaN and excluded from mean_fst")

    # Bin tables are keyed on an exact derived count, so they are only valid for
    # the sample size they were built at: 181 counts at n_hap=200, 155 at 172.
    # A mismatch would silently map every site to the wrong reference, so it is
    # a hard error here and again per sim in process_sim().
    man = os.path.join(bins_dir, "MANIFEST.txt")
    n_hap = None
    man_pops = None
    if os.path.exists(man):
        for line in open(man):
            if line.startswith("n_hap"):
                n_hap = int(line.split(":", 1)[1].strip())
            elif line.startswith("pop_ids"):
                man_pops = [t.strip() for t in
                            line.split(":", 1)[1].strip().split(",") if t.strip()]
    if n_hap is None:
        sys.exit(f"ERROR: {man} records no n_hap. Reference sets built before "
                 f"exact-count binning cannot be used; rebuild step 1.")
    if not man_pops:
        sys.exit(f"ERROR: {man} records no pop_ids. The population list cannot "
                 f"be inferred from the archive - see below.")
    log(f"reference set built at n_hap={n_hap}, pops {','.join(man_pops)}")
    # XP-EHH site set must match the null: the MANIFEST records whether the
    # bins keep target-monomorphic sites (k = 0 / n_hap). Mismatch = hard error.
    man_keep = "XPEHH_KEEP_TARGET_MONOMORPHIC : True" in open(man).read()
    if man_keep != C.XPEHH_KEEP_TARGET_MONOMORPHIC:
        sys.exit(f"ERROR: {man} XP-EHH bins keep target-monomorphic sites = "
                 f"{man_keep}, but cms_stats.XPEHH_KEEP_TARGET_MONOMORPHIC = "
                 f"{C.XPEHH_KEEP_TARGET_MONOMORPHIC}. Rebuild the bins "
                 f"(diag/rebuild_xpehh_bins_mono.py) or fix the setting.")
    log(f"XP-EHH keeps target-monomorphic sites: {man_keep}")

    work_dir = os.path.join(args.work_root, f"work_{args.model}_{batch_id}")
    if work_dir_holder is not None and not args.keep_scratch:
        work_dir_holder["path"] = work_dir
    for sub in ("dl", "tped", "rawstats", "scratch", "tsv"):
        os.makedirs(os.path.join(work_dir, sub), exist_ok=True)

    # ---- fetch ---------------------------------------------------------------
    if args.input_dir:
        src = args.input_dir
        log(f"reading tped zip + par_inputs metafile from {src}")
    else:
        src = os.path.join(work_dir, "dl")
        gcs = f"{GCS_BASE}/{args.model}"
        log("downloading tped zip + par_inputs metafile ...")
        sh(f"gcloud storage cp '{gcs}/*batch{batch_id}*tped*.zip' {src}/")
        sh(f"gcloud storage cp '{gcs}/*batch{batch_id}*par_inputs*.tsv' {src}/")

    zips = glob.glob(os.path.join(src, f"*batch{batch_id}*tped*.zip"))
    metas = glob.glob(os.path.join(src, f"*batch{batch_id}*par_inputs*.tsv"))
    if len(zips) != 1 or len(metas) != 1:
        sys.exit(f"expected 1 tped zip and 1 par_inputs tsv for batch {batch_id} "
                 f"in {src}, found {len(zips)} and {len(metas)}")
    # ---- TPED identity: (model, batch) from the file names -----------------
    # A TPED carries no pop or batch id in its content, so identity rests on
    # names; every one must agree with what was asked for.
    zb = os.path.basename(zips[0])
    zm = re.fullmatch(rf"{re.escape(args.model)}(_5mb)?_batch{batch_id}"
                      rf"_raw_tped(_\d{{4}}-\d{{2}}-\d{{2}})?\.zip", zb)
    if not zm:
        sys.exit(f"ERROR: IDENTITY: tped zip {zb} is not {args.model}[_5mb]_batch"
                 f"{batch_id}_raw_tped[_date].zip")
    if (zm.group(1) is not None) != (region_mb == 5):
        sys.exit(f"ERROR: IDENTITY: tped zip {zb} region tag disagrees with batch id "
                 f"region {region_mb} Mb")
    mb_ = os.path.basename(metas[0])
    if not re.match(rf"{re.escape(args.model)}(_5mb)?(\.par)?_batch{batch_id}_", mb_):
        sys.exit(f"ERROR: IDENTITY: metafile {mb_} does not belong to "
                 f"{args.model}:{batch_id}")
    meta = pd.read_csv(metas[0], sep="\t")
    log(f"metafile: {mb_} ({len(meta)} rows)")
    if "sim_batch_no" in meta.columns:
        other = set(meta["sim_batch_no"].astype(str)) - {str(batch_id)}
        if other:
            sys.exit(f"ERROR: IDENTITY: metafile {mb_} has sim_batch_no "
                     f"{sorted(other)[:5]} besides {batch_id}")

    # Index the archive WITHOUT extracting it. Workers extract their own
    # simulation and delete it again, so peak disk is bounded by
    # workers x one sim rather than by the whole batch.
    sim_groups, stray, dup = {}, [], []
    with zipfile.ZipFile(zips[0]) as zf:
        for member in zf.namelist():
            if not member.endswith(".tped"):
                continue
            m = TPED_MEMBER_RE.fullmatch(member)
            if not m:
                stray.append(member); continue
            g = sim_groups.setdefault(m.group(1), {})
            if m.group(2) in g:
                dup.append(member)
            g[m.group(2)] = member
    if stray or dup:
        sys.exit(f"ERROR: IDENTITY: tped zip {zb}: {len(stray)} member(s) not "
                 f"sel/sel.hap.<sim>_0_<pop>.tped {stray[:3]}, {len(dup)} duplicate "
                 f"(sim, pop) {dup[:3]}")
    if not sim_groups:
        sys.exit("no sel.hap.*.tped files in the archive")
    log(f"indexed {len(sim_groups)} sims in the archive (not extracted)")

    # The population list comes from the reference set, NOT from the archive.
    # cousins archives carry TPEDs for all six populations but only 5 and 6 have
    # nonzero Ne; pops 1-4 are written with 18,185 rows and ZERO haplotype
    # columns. Inferring pops from the archive would take those as outgroups,
    # attempt 5 XP-EHH pairs where only xpehh_5_vs_6_bin.csv exists, and fail
    # every sim with "position mismatch" - all 98 cousins batches.
    #
    # Taking the list from the MANIFEST also enforces the invariant that matters:
    # the batch is scored against a null built from the same populations.
    archive_pops = sorted(next(iter(sim_groups.values())))
    pop_ids = [p for p in man_pops if p in archive_pops]
    missing = [p for p in man_pops if p not in archive_pops]
    if missing:
        sys.exit(f"ERROR: reference set needs pops {man_pops} but the archive "
                 f"has only {archive_pops} (missing {missing})")
    dropped = [p for p in archive_pops if p not in man_pops]
    if dropped:
        log(f"ignoring {len(dropped)} population(s) present in the archive but "
            f"not in the reference set: {','.join(dropped)}")
    # Every sim must carry every reference pop. archive_pops above is read
    # from ONE sim; a sim missing an outgroup would otherwise get deldaf/FST
    # over fewer pops and max_xpehh over fewer pairs, silently.
    short = {sid: sorted(set(pop_ids) - set(pops))
             for sid, pops in sim_groups.items() if not set(pop_ids) <= set(pops)}
    if short:
        ex = ", ".join(f"{k}:{v}" for k, v in sorted(short.items(), key=lambda kv: int(kv[0]))[:10])
        sys.exit(f"ERROR: {len(short)} sim(s) lack reference pop TPEDs ({ex}); "
                 f"batch {args.model}:{batch_id} NOT written")
    sim_groups = {sid: {p: f for p, f in pops.items() if p in pop_ids}
                  for sid, pops in sim_groups.items()}
    alt_pops = [p for p in pop_ids if p != target_pop]
    pairs = [f"{target_pop}_vs_{p}" for p in alt_pops]
    log(f"{len(sim_groups)} sims | pops {pop_ids} | pairs {pairs}")

    # Fail before computing anything: one TPED is enough to know the sample size.
    with zipfile.ZipFile(zips[0]) as zf:
        _pm = next(iter(sim_groups.values()))[target_pop]
        _probe_f = os.path.join(work_dir, "probe.tped")
        with zf.open(_pm) as _s, open(_probe_f, "wb") as _o:
            shutil.copyfileobj(_s, _o)
    probe = C.parse_tped(_probe_f)[3]
    os.remove(_probe_f)
    if probe != n_hap:
        sys.exit(f"ERROR: batch has {probe} haplotypes per population but the "
                 f"reference set at {bins_dir} was built at n_hap={n_hap}. "
                 f"Per-count bins are not transferable between sample sizes.")

    bins = load_bins(bins_dir, pairs, flavours)

    # ---- raw stats in parallel ------------------------------------------------
    # Two sims in ~1,028,000 have no causal row at all and are dropped by
    # decision (2026-10-15). Skipping them explicitly and logging it keeps the
    # count reproducible instead of silently differing.
    drop = set(DROP_SIMS.get((args.model, batch_id), set()))
    if args.drop_list:
        dl = pd.read_csv(args.drop_list, sep="\t", dtype={"sim_batch_no": str})
        need = {"model", "sim_batch_no", "sim_id"}
        if not need <= set(dl.columns):
            sys.exit(f"ERROR: --drop-list needs columns {sorted(need)}")
        # keyed on (model, batch): batch ids repeat across models
        mine = dl[(dl.model == args.model) & (dl.sim_batch_no == str(batch_id))]
        drop |= set(mine.sim_id.astype(int).astype(str))
        log(f"drop list {os.path.basename(args.drop_list)}: {len(mine)} sim(s) "
            f"listed for ({args.model}, {batch_id})")
    if drop:
        present = sorted(drop & set(sim_groups), key=int)
        for sid in present:
            del sim_groups[sid]
        log(f"dropping {len(present)} sim(s) (no causal row, or causal DAF "
            f"outside 0.05-0.95 per --drop-list): {','.join(present)}")
        absent = sorted(drop - set(present), key=int)
        if absent:
            log(f"WARNING: {len(absent)} listed sim(s) not in the archive: "
                f"{','.join(absent[:20])}")

    # Metafile causal DAF per sim, for check_target_label. Archived metafiles
    # hold selpos rows; the replacement par_inputs one row per sim. If the
    # metafile has a daf column, every computed sim must have exactly one row.
    expected = {}
    if "daf" in meta.columns:
        mc = meta[meta["selpos"] == 1] if "selpos" in meta.columns else meta
        cnt = mc["sim_id"].astype(int).astype(str).value_counts()
        dup = sorted(set(cnt[cnt > 1].index) & set(sim_groups), key=int)
        miss = sorted(set(sim_groups) - set(cnt.index), key=int)
        if dup or miss:
            sys.exit(f"ERROR: {LABEL_FAIL}: metafile causal rows: "
                     f"{len(dup)} sim(s) with >1 ({','.join(dup[:10])}), "
                     f"{len(miss)} with none ({','.join(miss[:10])})")
        expected = dict(zip(mc["sim_id"].astype(int).astype(str),
                            mc["daf"].astype(float)))
    else:
        log("WARNING: metafile has no daf column; label check uses the causal "
            "MAF gate only")
    tasks = [(sid, pops, target_pop, work_dir, flavours, n_hap, zips[0],
              selpos, expected.get(sid), tuple(pop_ids))
             for sid, pops in sorted(sim_groups.items(), key=lambda kv: int(kv[0]))]
    log(f"computing raw statistics with {args.workers} workers ...")
    ok, failed = 0, []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for fut in as_completed([ex.submit(process_sim, t) for t in tasks]):
            good, sid, msg = fut.result()
            if good:
                ok += 1
                if ok % 50 == 0:
                    log(f"  {ok}/{len(tasks)} sims")
            else:
                failed.append((sid, msg))
    log(f"raw stage: {ok} ok, {len(failed)} failed")
    # Workers read the archive directly, so it can only go once they finish.
    if not args.input_dir:
        try:
            os.remove(zips[0])   # never delete the user's local input copy
        except OSError:
            pass
    for sid, msg in failed[:10]:
        log(f"  FAILED sim {sid}: {msg}")
    guard_bad = sorted((x for x in failed if x[1].startswith(GUARD_FAIL)),
                       key=lambda x: int(x[0]))
    if guard_bad:
        for sid, msg in guard_bad[:20]:
            log(f"  {sid}: {msg}")
        sys.exit(f"ERROR: {len(guard_bad)} of {len(tasks)} sim(s) fail an input "
                 f"guard; batch {args.model}:{batch_id} NOT written")
    label_bad = sorted((x for x in failed if x[1].startswith(LABEL_FAIL)),
                       key=lambda x: int(x[0]))
    if label_bad:
        for sid, msg in label_bad[:20]:
            log(f"  {sid}: {msg}")
        sys.exit(f"ERROR: {len(label_bad)} of {len(tasks)} sim(s) fail the target "
                 f"label check; batch {args.model}:{batch_id} NOT written")
    if failed:
        sys.exit(f"ERROR: {len(failed)} of {len(tasks)} sim(s) failed the raw "
                 f"stage; batch {args.model}:{batch_id} NOT written")

    # ---- normalize, collate ---------------------------------------------------
    log("normalizing and collating ...")
    meta_idx = (meta.set_index(meta["sim_id"].astype(str))
                if "sim_id" in meta.columns else pd.DataFrame())

    # Normalization runs in the same process pool as the raw stage: run
    # serially it took 2.5 of 12.9 min on a gravel 3 Mb batch with the other
    # cores idle (2026-10-16). Per-sim logic is unchanged (finalize_sim).
    meta_rows = {sid: meta_idx.loc[sid].to_dict() for sid in sim_groups
                 if sid in meta_idx.index}
    ftasks = [(sid, work_dir, flavours, bins, n_hap, pairs, batch_id, selpos,
               args.model, meta_rows.get(sid))
              for sid in sorted(sim_groups, key=int)]
    written = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for fut in as_completed([ex.submit(finalize_sim, t) for t in ftasks]):
            ok_, msgs = fut.result()
            written += ok_
            for m_ in msgs:
                log(m_)

    log(f"wrote {written} per-sim TSVs")
    if written != len(tasks):
        sys.exit(f"ERROR: wrote {written} of {len(tasks)} sims; "
                 f"batch {args.model}:{batch_id} NOT written")

    # ---- quick self-check -----------------------------------------------------
    verify(work_dir, pairs, selpos, flavours)

    # ---- zip -------------------------------------------------------------------
    os.makedirs(out_dir, exist_ok=True)
    final = os.path.join(out_dir,
                         f"{args.model}_batch{batch_id}_cms_stats_recomputed.zip")
    tmp = final + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for f in sorted(glob.glob(os.path.join(work_dir, "tsv", "*.tsv"))):
            zf.write(f, os.path.basename(f))
        for f in sorted(glob.glob(os.path.join(bins_dir, "*"))):
            zf.write(f, os.path.join("bins", os.path.basename(f)))
    os.replace(tmp, final)
    log(f"archive: {final}")

    if not args.keep_scratch:
        shutil.rmtree(work_dir, ignore_errors=True)
    log("done")


# ==============================================================================
def verify(work_dir, pairs, selpos, flavours):
    """Checks on the things most likely to be silently wrong."""
    log("--- verification ---")
    tsvs = sorted(glob.glob(os.path.join(work_dir, "tsv", "*.tsv")))
    if not tsvs:
        log("  nothing to verify")
        return

    sample = tsvs[min(1, len(tsvs) - 1)]      # not sim 0: see CLAUDE.md
    sid = re.search(r"_cms_stats_(\d+)\.tsv$", sample).group(1)

    # per-pair raw row counts: unequal counts are exactly the condition that
    # used to break row-order alignment
    for pair in pairs:
        f = os.path.join(work_dir, "rawstats", f"{sid}_{pair}.xpehh.out")
        if os.path.exists(f):
            n = sum(1 for _ in open(f)) - 1
            log(f"  sim {sid} raw {pair}: {n:,} scored sites")

    allrows = pd.concat([pd.read_csv(f, sep="\t") for f in tsvs[:50]],
                        ignore_index=True)
    for col in out_cols(flavours):
        if col in ("sim_batch_no", "sim_id", "pos", "selpos",
                   "deri_gen", "sel_gen", "s"):
            continue
        log(f"  {col:<22} missing {100 * allrows[col].isna().mean():5.2f}%")

    sel = allrows[allrows["pos"] == selpos]
    if not sel.empty:
        log(f"  selected sites sampled: {len(sel)}, "
            f"daf {sel['daf'].min():.3f}-{sel['daf'].max():.3f}")
        for flav in flavours:
            n = sel[f"norm_ihs_{flav}"].isna().sum()
            log(f"  selected sites with no iHS ({flav}): {n}/{len(sel)}")


if __name__ == "__main__":
    main()


# NOTE on cleanup: main() removes work_dir on the success path. The wrapper
# below guarantees it on EVERY path, including a crash or SIGTERM from a
# walltime kill. Without it a failed task leaves the whole extracted archive
# behind - up to 69 GB for medina-munoz 5 Mb - and 32 concurrent tasks would
# fill node-local disk within a few failures.

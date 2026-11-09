#!/usr/bin/env python3
"""Simulate and score one neutral batch (one row of the neutral plan TSV).

Per sim:
  1. cosi2 on the reference set's own resolved.par (bins/<combo>/resolved.par,
     the exact par the null was built from), no sweep line, production flags.
     Seed = seed0 + sim*16 + attempt. Any failed attempt (timeout, non-zero
     exit, missing TPED) is retried with the next seed; every attempt is logged.
  2. flip 0<->1 (derived = 1), as every pipeline does.
  3. score with 02_recompute_batch.process_sim + finalize_sim, imported from
     the production file, so the neutral statistics are computed by the same
     code as the 1,010,000-sim corpus. Two hooks are replaced, both because a
     neutral sim has no causal variant:
       check_target_label -> skipped (it requires a causal site at selpos)
       out_cols           -> the neutral schema, without the sweep columns
                             deri_gen, sel_gen, s, selpos (raw_sim_neut_2026)
  4. TPEDs are deleted. Only the seed is kept: the sim is regenerated from it.

Outputs (out/<model>/):
  <par_base>.par_batch<id>_cms_stats_all_<date>.zip    one TSV per sim
  <par_base>.par_batch<id>_neutral_seeds_<date>.tsv    sim_id, seed, attempts
and state/<model>__<batch>.json, written last, only if every sim succeeded.

Usage: run_neutral_batch.py --batches PLAN.tsv --row IDX --out-root DIR
                            --state-root DIR [--n-sims N] [--workers W]
"""
import argparse, csv, glob, hashlib, importlib.util, json, os, shutil, subprocess
import sys, time, zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

R = os.path.dirname(os.path.abspath(__file__))   # pipeline/: cms_stats, 02
sys.path.insert(0, R)
import cms_stats as C  # noqa: E402

_spec = importlib.util.spec_from_file_location("r02", f"{R}/02_recompute_batch.py")
R02 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R02)

NEUT_COLS = ["sim_batch_no", "sim_id", "pos", "mean_fst", "deldaf", "daf", "maf",
             "norm_nsl", "norm_ihh12", "norm_ihs_binom", "norm_delihh_binom",
             "max_xpehh"]
R02.check_target_label = lambda *a, **k: None
R02.out_cols = lambda flavours: NEUT_COLS

ATTEMPT_TIMEOUT, SIM_BUDGET, MAX_ATTEMPTS = 120.0, 600.0, 16
COSI_FLAGS = ["-v", "-g", "--genmapRandomRegions", "--nTriesPerSim", "5000",
              "--drop-singletons", ".25", "-n", "1", "-M"]


def log(m):
    print(f"[{datetime.now():%F %T}] {m}", flush=True)


def simulate(par, pop_ids, sim_dir, seed0, sid):
    """Run cosi2 until one attempt yields every pop's TPED; flip into
    sim_dir/sel/sel.hap.<sid>_0_<pop>.tped. Returns (seed, attempts, log)."""
    prefix = os.path.join(sim_dir, "raw")
    env = dict(os.environ, COSI_NEWSIM="1", COSI_MAXATTEMPTS="5000")
    t0, tries = time.time(), []
    for a in range(MAX_ATTEMPTS):
        spent = time.time() - t0
        if spent >= SIM_BUDGET:
            break
        seed = seed0 + sid * 16 + a
        for f in glob.glob(f"{prefix}_*"):
            os.remove(f)
        a0 = time.time()
        try:
            subprocess.run(["coalescent", "-p", par, *COSI_FLAGS, "--tped", prefix,
                            "-r", str(seed)], check=True, capture_output=True,
                           env=env, cwd=sim_dir, timeout=min(ATTEMPT_TIMEOUT, SIM_BUDGET - spent))
            outcome = "ok"
        except subprocess.TimeoutExpired:
            outcome = "timeout"
        except subprocess.CalledProcessError as e:
            outcome = f"exit{e.returncode}"
        if outcome == "ok" and not all(os.path.exists(f"{prefix}_0_{p}.tped") for p in pop_ids):
            outcome = "no_tped"
        tries.append((seed, outcome, round(time.time() - a0, 1)))
        if outcome == "ok":
            os.makedirs(os.path.join(sim_dir, "sel"), exist_ok=True)
            for p in pop_ids:
                C.flip_tped(f"{prefix}_0_{p}.tped",
                            os.path.join(sim_dir, "sel", f"sel.hap.{sid}_0_{p}.tped"))
            for f in glob.glob(f"{prefix}_*"):
                os.remove(f)
            return seed, a + 1, tries
    raise RuntimeError(f"sim {sid}: no successful cosi2 attempt: {tries}")


def score(sid, tped_dir, work_dir, target_pop, pop_ids, flavours, n_hap, bins,
          pairs, batch_id, selpos, member_prefix):
    """Score flipped TPEDs with 02's own code. TPEDs are packed into a STORED
    zip because process_sim reads (and identity-checks) zip members."""
    zpath = os.path.join(tped_dir, "sim.zip")
    members = {}
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_STORED) as zf:
        for p in pop_ids:
            m = f"sel/sel.hap.{sid}_0_{p}.tped"
            f = os.path.join(tped_dir, m)
            zf.write(f, m); os.remove(f)
            members[str(p)] = m
    ok, _, msg = R02.process_sim((str(sid), members, str(target_pop), work_dir,
                                  flavours, n_hap, zpath, selpos, None,
                                  tuple(str(p) for p in pop_ids)))
    os.remove(zpath)
    if not ok:
        raise RuntimeError(f"sim {sid}: process_sim: {msg}")
    wrote, msgs = R02.finalize_sim((str(sid), work_dir, flavours, bins, n_hap, pairs,
                                    batch_id, selpos, member_prefix, None))
    for f in glob.glob(os.path.join(work_dir, "rawstats", f"{sid}.raw.tsv")) + \
            glob.glob(os.path.join(work_dir, "rawstats", f"{sid}_*.xpehh.out")):
        os.remove(f)
    if wrote != 1:
        raise RuntimeError(f"sim {sid}: finalize_sim: {msgs}")


def one_sim(task):
    (sid, seed0, par, pop_ids, target_pop, work_dir, flavours, n_hap, bins, pairs,
     batch_id, selpos, member_prefix) = task
    sim_dir = os.path.join(work_dir, "sims", str(sid))
    os.makedirs(sim_dir, exist_ok=True)
    t0 = time.time()
    try:
        seed, attempts, tries = simulate(par, pop_ids, sim_dir, seed0, sid)
        t1 = time.time()
        score(sid, sim_dir, work_dir, target_pop, pop_ids, flavours, n_hap, bins,
              pairs, batch_id, selpos, member_prefix)
        return dict(sim_id=sid, seed=seed, attempts=attempts,
                    failed_attempts=";".join(f"{s}:{o}" for s, o, _ in tries if o != "ok"),
                    sim_secs=round(t1 - t0, 1), score_secs=round(time.time() - t1, 1))
    finally:
        shutil.rmtree(sim_dir, ignore_errors=True)


def read_manifest(bins_dir):
    txt = open(os.path.join(bins_dir, "MANIFEST.txt")).read()
    n_hap = pops = None
    for line in txt.splitlines():
        if line.startswith("n_hap"):
            n_hap = int(line.split(":", 1)[1])
        elif line.startswith("pop_ids"):
            pops = [int(x) for x in line.split(":", 1)[1].split(",")]
    keep = "XPEHH_KEEP_TARGET_MONOMORPHIC : True" in txt
    assert n_hap and pops, bins_dir
    assert keep == C.XPEHH_KEEP_TARGET_MONOMORPHIC, f"{bins_dir}: XP-EHH bin site set mismatch"
    return n_hap, pops, hashlib.md5(txt.encode()).hexdigest()


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batches", required=True, help="neutral plan TSV (workflow/plan.py)")
    ap.add_argument("--row", type=int, required=True, help="plan idx")
    ap.add_argument("--n-sims", type=int, help="override (tests only)")
    ap.add_argument("--first-sim", type=int, default=0, help="tests only")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--state-root", required=True)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    a = ap.parse_args()

    row = [r for r in csv.DictReader(open(a.batches), delimiter="\t")
           if int(r["idx"]) == a.row][0]
    model, batch, target_pop = row["model"], row["batch"], int(row["target_pop"])
    state = os.path.join(a.state_root, f"{model}__{batch}.json")
    if os.path.exists(state):
        log(f"{model}:{batch} already done ({state})"); return
    C.check_tools(["coalescent", "selscan", "ihsbin", "xpehhbin"])
    flavours = C.ihs_flavours(False)
    assert list(flavours) == ["binom"], flavours

    tp, mb, selpos = C.decode_batch_id(batch)
    assert (int(tp), mb) == (target_pop, int(row["region_mb"])), (batch, tp, mb)
    n_hap, pops, man_md5 = read_manifest(row["bins_dir"])
    assert pops == [int(x) for x in row["pop_ids"].split(",")], (pops, row["pop_ids"])
    par = os.path.join(row["bins_dir"], "resolved.par")
    length = [int(l.split()[1]) for l in open(par) if l.split()[:1] == ["length"]]
    assert length == [mb * 1_000_000], (par, length)
    pairs = [f"{target_pop}_vs_{p}" for p in pops if p != target_pop]
    bins = R02.load_bins(row["bins_dir"], pairs, flavours)

    n = a.n_sims or int(row["n_sims"])
    sims = list(range(a.first_sim, a.first_sim + n))
    assert sims[-1] * 16 + 15 < 20_000, "seed block overflow"
    member_prefix = f"{row['par_base']}.par"
    log(f"{model}:{batch} | {row['combo']} | {n} sims | seed0 {row['seed0']} | "
        f"n_hap {n_hap} | pairs {pairs} | {a.workers} workers")

    scratch = os.environ.get("TMPDIR") or "/tmp"
    work_dir = os.path.join(scratch, f"neut_{model}_{batch}_{os.getpid()}")
    for sub in ("rawstats", "scratch", "tsv", "sims"):
        os.makedirs(os.path.join(work_dir, sub), exist_ok=True)
    t0 = time.time()
    try:
        tasks = [(s, int(row["seed0"]), par, pops, target_pop, work_dir, flavours,
                  n_hap, bins, pairs, batch, selpos, member_prefix) for s in sims]
        recs, fails = [], []
        with ProcessPoolExecutor(a.workers) as ex:
            futs = {ex.submit(one_sim, t): t[0] for t in tasks}
            for f in as_completed(futs):
                try:
                    recs.append(f.result())
                    if len(recs) % 100 == 0:
                        log(f"  {len(recs)}/{n} sims ({time.time() - t0:.0f} s)")
                except Exception as e:
                    fails.append((futs[f], repr(e)[:400]))
        for s, e in fails[:10]:
            log(f"  FAILED sim {s}: {e}")
        if fails:
            sys.exit(f"ERROR: {len(fails)} of {n} sims failed; batch NOT written")

        recs.sort(key=lambda r: r["sim_id"])
        tsvs = sorted(glob.glob(os.path.join(work_dir, "tsv", "*.tsv")))
        want = {f"{member_prefix}_batch{batch}_cms_stats_{s}.tsv" for s in sims}
        assert {os.path.basename(t) for t in tsvs} == want, "TSV set != sims"
        for t in tsvs:
            with open(t) as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                assert hdr == NEUT_COLS, (t, hdr)
                assert fh.readline(), f"{t}: no sites"
        assert len({r["seed"] for r in recs}) == n, "duplicate seeds"

        out = os.path.join(a.out_root, model)
        os.makedirs(out, exist_ok=True)
        zfinal = os.path.join(out, row["zip_name"])
        with zipfile.ZipFile(zfinal + ".tmp", "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            for t in tsvs:
                zf.write(t, os.path.basename(t))
        with zipfile.ZipFile(zfinal + ".tmp") as zf:
            assert zf.testzip() is None and len(zf.namelist()) == n
        os.replace(zfinal + ".tmp", zfinal)
        sfinal = os.path.join(out, row["seeds_name"])
        with open(sfinal + ".tmp", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sim_batch_no", "sim_id", "seed", "attempts",
                                               "failed_attempts", "sim_secs", "score_secs"],
                               delimiter="\t")
            w.writeheader()
            for r in recs:
                w.writerow(dict(sim_batch_no=batch, **r))
        os.replace(sfinal + ".tmp", sfinal)

        os.makedirs(a.state_root, exist_ok=True)
        info = dict(row, n_written=n, first_sim=a.first_sim, zip=zfinal, seeds=sfinal,
                    zip_bytes=os.path.getsize(zfinal), n_hap=n_hap,
                    bins_manifest_md5=man_md5, resolved_par_md5=md5(par),
                    code_md5={f: md5(os.path.join(R, f)) for f in
                              ("cms_stats.py", "02_recompute_batch.py",
                               "run_neutral_batch.py")},
                    tools=C.tool_versions(), cosi_flags=" ".join(COSI_FLAGS),
                    reseeded=sum(r["attempts"] > 1 for r in recs),
                    secs=round(time.time() - t0, 1), job=os.environ.get("SLURM_JOB_ID"),
                    finished=datetime.now().isoformat(timespec="seconds"))
        with open(state + ".tmp", "w") as fh:
            json.dump(info, fh, indent=1, default=str)
        os.replace(state + ".tmp", state)
        log(f"done: {n} sims, {info['reseeded']} reseeded, {info['secs']} s -> {zfinal}")
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    main()

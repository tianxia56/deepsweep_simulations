#!/usr/bin/env python3
"""Simulate one batch of selected-sweep replicates with cosi2.

Replaces the production stage 02cs_sel.sh + finish_batch.py with the same
rules, so a batch with the same seed0 regenerates the same TPEDs:

  par      base model .par (recomb_file made absolute) + the sweep line
  cosi2    COSI_NEWSIM=1 COSI_MAXATTEMPTS=5000 coalescent -p <par> -v -g
           --genmapRandomRegions --nTriesPerSim 5000 --drop-singletons .25
           --tped <prefix> -n 1 -M -r <seed>
  seed     seed0 + sim_id*16 + attempt, attempt = 1..15 (stride 16)
  limits   120 s per attempt, 600 s per sim; a sim that exhausts either is
           SKIPPED and logged, never retried with the same seed
  replay   --replay TSV (sim_id, seed): run exactly those sims, each once at
           its recorded seed and with no time limit, so a slower machine
           cannot change which attempt succeeds. Used to regenerate the
           2026oct corpus (provenance/sweep_sims_kept_2026oct.tsv).
  coding   flipped 0<->1, so derived = 1 (cms_stats.flip_tped)
  params   deri_gen, s, sel_gen from the LAST sweep_mult_standing line cosi2
           wrote to COSI_SAVE_SAMPLED for the successful attempt (fields 3, 4,
           8, as 08_output_tped.sh parsed them)

Outputs in --out-dir, named as 02_recompute_batch.py's identity check expects:

  <par_base>_batch<id>_raw_tped.zip          sel/sel.hap.<sim>_0_<pop>.tped
                                             + sim_log/seeds.tsv
  <par_base>.par_batch<id>_par_inputs.tsv    one row per simulated sim:
        sim_batch_no sim_id pos daf n_hap deri_gen sel_gen s seed attempt selpos
        (daf = causal DAF in the target pop, NA if no site at selpos)
  <par_base>_batch<id>_drop.tsv              sims to leave out of scoring:
        no causal row, or causal DAF outside [0.05, 0.95] (the causal MAF
        gate of 02's check_target_label, same rounding)
  <par_base>_batch<id>_summary.json          counts; written LAST

Usage:
  simulate_sweep_batch.py --model gravel --par models/gravel.par --region-mb 3
      --target-pop 1 --pop-ids 1,2,3 --batch-id 1301000000 --n-sims 100
      --seed0 2147483648 --sweep 'pop_event sweep_mult_standing ...'
      --out-dir DIR [--workers N]
"""
import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cms_stats as C  # noqa: E402

ATTEMPT_TIMEOUT, SIM_BUDGET, MAX_ATTEMPTS, SEED_STRIDE = 120.0, 600.0, 15, 16
COSI_FLAGS = ["-v", "-g", "--genmapRandomRegions", "--nTriesPerSim", "5000",
              "--drop-singletons", ".25", "-n", "1", "-M"]
DAF_MIN, DAF_MAX = 0.05, 0.95


def log(m):
    print(f"[{datetime.now():%F %T}] {m}", flush=True)


def write_par(base_par, sweep, region_mb, out):
    """Base .par with recomb_file absolute, then the sweep line appended."""
    par_dir = os.path.dirname(os.path.abspath(base_par))
    lines, length = [], None
    for line in open(base_par):
        tok = line.split()
        if tok and tok[0] == "length":
            length = int(tok[1])
        if tok and tok[0] == "recomb_file" and not os.path.isabs(tok[1]):
            recom = os.path.join(par_dir, tok[1])
            if not os.path.exists(recom):
                sys.exit(f"recomb_file not found: {recom}")
            line = f"recomb_file {recom}\n"
        lines.append(line)
    if length != region_mb * 1_000_000:
        sys.exit(f"{base_par} has length {length}, but --region-mb is {region_mb}")
    if not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    with open(out, "w") as fh:
        fh.writelines(lines)
        fh.write(sweep.rstrip("\n") + "\n")
    return out


def sampled_params(path):
    last = None
    for line in open(path):
        p = line.split()
        if p[:1] == ["sweep_mult_standing"]:
            last = p
    if last is None:
        raise RuntimeError(f"{path}: no sweep_mult_standing line")
    return dict(deri_gen=round(float(last[3])), s=round(float(last[4]), 6),
                sel_gen=round(float(last[8])))


def causal_daf(tped, selpos):
    for line in open(tped):
        p = line.split()
        if int(p[3]) == selpos:
            g = p[4:]
            assert set(g) <= {"0", "1"}, f"{tped}: non-0/1 allele at selpos"
            return g.count("1"), len(g)
    return None, None


def one_sim(task):
    """-> dict(sim_id, ok, seed, attempt, params, tries, tpeds). TPEDs are left
    flipped in sim_dir/sel/ for the parent to pack."""
    sid, par, pop_ids, sim_dir, seed0, fixed_seed = task
    os.makedirs(sim_dir, exist_ok=True)
    prefix = os.path.join(sim_dir, "raw")
    sampled = os.path.join(sim_dir, "sampled.txt")
    t0, tries = time.time(), []
    attempts = ([fixed_seed - seed0 - sid * SEED_STRIDE] if fixed_seed is not None
                else range(1, MAX_ATTEMPTS + 1))
    for attempt in attempts:
        spent = time.time() - t0
        if fixed_seed is None and spent >= SIM_BUDGET:
            break
        seed = seed0 + sid * SEED_STRIDE + attempt
        for f in glob.glob(f"{prefix}_*") + [sampled]:
            if os.path.exists(f):
                os.remove(f)
        env = dict(os.environ, COSI_NEWSIM="1", COSI_MAXATTEMPTS="5000",
                   COSI_SAVE_SAMPLED=sampled)
        a0 = time.time()
        try:
            subprocess.run(["coalescent", "-p", par, *COSI_FLAGS, "--tped", prefix,
                            "-r", str(seed)], check=True, capture_output=True,
                           env=env, cwd=sim_dir,
                           timeout=(None if fixed_seed is not None
                                    else min(ATTEMPT_TIMEOUT, SIM_BUDGET - spent)))
            outcome = "ok"
        except subprocess.TimeoutExpired:
            outcome = "timeout"
        except subprocess.CalledProcessError as e:
            outcome = f"exit{e.returncode}"
        written = sorted(glob.glob(f"{prefix}_0_*.tped"))
        if outcome == "ok" and not all(os.path.exists(f"{prefix}_0_{p}.tped")
                                       for p in pop_ids):
            outcome = "no_tped"
        tries.append((attempt, seed, outcome, round(time.time() - a0, 1)))
        if outcome == "ok":
            params = sampled_params(sampled)
            os.makedirs(os.path.join(sim_dir, "sel"), exist_ok=True)
            tpeds = {}
            for f in written:
                pop = os.path.basename(f)[len("raw_0_"):-len(".tped")]
                m = f"sel/sel.hap.{sid}_0_{pop}.tped"
                C.flip_tped(f, os.path.join(sim_dir, m))
                os.remove(f)
                tpeds[pop] = m
            return dict(sim_id=sid, ok=True, seed=seed, attempt=attempt,
                        params=params, tries=tries, tpeds=tpeds)
    return dict(sim_id=sid, ok=False, tries=tries)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--par", required=True, help="base model .par, no sweep line")
    ap.add_argument("--region-mb", type=int, required=True)
    ap.add_argument("--target-pop", type=int, required=True)
    ap.add_argument("--pop-ids", required=True)
    ap.add_argument("--batch-id", required=True)
    ap.add_argument("--n-sims", type=int, required=True)
    ap.add_argument("--seed0", type=int, required=True)
    ap.add_argument("--sweep", required=True, help="cosi2 sweep line")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--first-sim", type=int, default=0, help="tests only")
    ap.add_argument("--replay", default=None,
                    help="TSV with sim_id, seed: simulate exactly these sims at these seeds")
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    a = ap.parse_args()

    pop_ids = [int(p) for p in a.pop_ids.split(",")]
    tp, mb, selpos = C.decode_batch_id(a.batch_id)
    if (int(tp), mb) != (a.target_pop, a.region_mb):
        sys.exit(f"batch id {a.batch_id} encodes pop {tp}, {mb} Mb; "
                 f"asked pop {a.target_pop}, {a.region_mb} Mb")
    if a.target_pop not in pop_ids:
        sys.exit(f"target pop {a.target_pop} not in {pop_ids}")
    if f" {a.target_pop} " not in f" {a.sweep} ":
        sys.exit(f"sweep line does not name target pop {a.target_pop}: {a.sweep}")
    sims = list(range(a.first_sim, a.first_sim + a.n_sims))
    fixed = {}
    if a.replay:
        with open(a.replay) as fh:
            fixed = {int(r["sim_id"]): int(r["seed"]) for r in csv.DictReader(fh, delimiter="\t")}
        sims = sorted(fixed)
        for sid, sd in fixed.items():
            if not 1 <= sd - a.seed0 - sid * SEED_STRIDE <= MAX_ATTEMPTS:
                sys.exit(f"replay sim {sid}: seed {sd} is not seed0 + sim*16 + attempt")
    if a.seed0 < 1 or a.seed0 + (sims[-1] + 1) * SEED_STRIDE >= 2**32:
        sys.exit(f"seed0 {a.seed0}: seeds leave cosi2's 32-bit range")
    C.check_tools(["coalescent"])

    par_base = os.path.basename(a.par)[:-len(".par")]
    stem = f"{par_base}_batch{a.batch_id}"
    os.makedirs(a.out_dir, exist_ok=True)
    zfinal = os.path.join(a.out_dir, f"{stem}_raw_tped.zip")
    summ_f = os.path.join(a.out_dir, f"{stem}_summary.json")
    if os.path.exists(summ_f):
        log(f"already done: {summ_f}")
        return
    scratch = os.environ.get("TMPDIR") or "/tmp"
    work = os.path.join(scratch, f"sweep_{a.model}_{a.batch_id}_{os.getpid()}")
    os.makedirs(work)
    par = write_par(a.par, a.sweep, a.region_mb, os.path.join(work, "sweep.par"))
    log(f"{a.model}:{a.batch_id} | pop {a.target_pop} of {pop_ids} | {a.region_mb} Mb "
        f"| selpos {selpos:,} | {a.n_sims} sims | seed0 {a.seed0} | {a.workers} workers")
    log(f"sweep: {a.sweep}")

    rows, drops, seed_log, skipped = [], [], [], []
    t0 = time.time()
    try:
        with zipfile.ZipFile(zfinal + ".tmp", "w", zipfile.ZIP_DEFLATED,
                             compresslevel=9) as zf, \
                ProcessPoolExecutor(a.workers) as ex:
            futs = [ex.submit(one_sim, (s, par, pop_ids,
                                        os.path.join(work, "sims", str(s)), a.seed0,
                                        fixed.get(s)))
                    for s in sims]
            for n_done, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                sid = r["sim_id"]
                sim_dir = os.path.join(work, "sims", str(sid))
                seed_log += [dict(sim_id=sid, attempt=at, seed=sd, outcome=o, secs=t)
                             for at, sd, o, t in r["tries"]]
                if not r["ok"]:
                    if fixed:
                        sys.exit(f"ERROR: replay sim {sid} failed at its recorded seed: {r['tries']}")
                    skipped.append(sid)
                    shutil.rmtree(sim_dir, ignore_errors=True)
                    continue
                k, n = causal_daf(os.path.join(sim_dir, r["tpeds"][str(a.target_pop)]),
                                  selpos)
                for m in sorted(r["tpeds"].values()):
                    zf.write(os.path.join(sim_dir, m), m)
                shutil.rmtree(sim_dir, ignore_errors=True)
                daf = None if k is None else k / n
                rows.append(dict(sim_batch_no=a.batch_id, sim_id=sid, pos=selpos,
                                 daf=daf, n_hap=n, **r["params"], seed=r["seed"],
                                 attempt=r["attempt"], selpos=1))
                reason = ("no_causal_row" if daf is None else
                          f"daf_below_{DAF_MIN}" if round(daf, 6) < DAF_MIN else
                          f"daf_above_{DAF_MAX}" if round(daf, 6) > DAF_MAX else None)
                if reason:
                    drops.append(dict(model=a.model, sim_batch_no=a.batch_id,
                                      sim_id=sid, daf=daf, reason=reason))
                if n_done % 25 == 0:
                    log(f"  {n_done}/{len(sims)} done ({len(skipped)} skipped, "
                        f"{time.time() - t0:.0f} s)")
            seed_log.sort(key=lambda d: (d["sim_id"], d["attempt"]))
            buf = "sim_id\tattempt\tseed\toutcome\tsecs\n" + "".join(
                f"{d['sim_id']}\t{d['attempt']}\t{d['seed']}\t{d['outcome']}\t{d['secs']}\n"
                for d in seed_log)
            zf.writestr("sim_log/seeds.tsv", buf)
        if not rows:
            sys.exit("ERROR: no sim succeeded; nothing written")
        with zipfile.ZipFile(zfinal + ".tmp") as zf:
            assert zf.testzip() is None
            got = {int(m.split(".")[2].split("_")[0]) for m in zf.namelist()
                   if m.startswith("sel/")}
        assert got == {r["sim_id"] for r in rows}, "zip sims != simulated sims"
        os.replace(zfinal + ".tmp", zfinal)

        rows.sort(key=lambda d: d["sim_id"])
        meta_f = os.path.join(a.out_dir, f"{par_base}.par_batch{a.batch_id}_par_inputs.tsv")
        cols = ["sim_batch_no", "sim_id", "pos", "daf", "n_hap", "deri_gen",
                "sel_gen", "s", "seed", "attempt", "selpos"]
        with open(meta_f, "w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({**r, "daf": "NA" if r["daf"] is None else repr(r["daf"])})
        drop_f = os.path.join(a.out_dir, f"{stem}_drop.tsv")
        with open(drop_f, "w", newline="") as fh:
            w = csv.DictWriter(fh, ["model", "sim_batch_no", "sim_id", "daf", "reason"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(sorted(drops, key=lambda d: d["sim_id"]))
        summ = dict(model=a.model, batch_id=a.batch_id, par_base=par_base,
                    target_pop=a.target_pop, pop_ids=pop_ids, region_mb=a.region_mb,
                    sweep=a.sweep, seed0=a.seed0, asked=len(sims), replay=bool(fixed),
                    simulated=len(rows), skipped=sorted(skipped),
                    dropped_for_scoring=len(drops), to_score=len(rows) - len(drops),
                    attempts=len(seed_log),
                    timeouts=sum(d["outcome"] == "timeout" for d in seed_log),
                    cosi_flags=" ".join(COSI_FLAGS), zip=zfinal, par_inputs=meta_f,
                    drop_list=drop_f, secs=round(time.time() - t0, 1),
                    job=os.environ.get("SLURM_JOB_ID"),
                    finished=datetime.now().isoformat(timespec="seconds"))
        with open(summ_f + ".tmp", "w") as fh:
            json.dump(summ, fh, indent=1)
        os.replace(summ_f + ".tmp", summ_f)
        log(f"done: {len(rows)} simulated, {len(skipped)} skipped, {len(drops)} "
            f"outside the causal DAF gate -> {zfinal}")
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if os.path.exists(zfinal + ".tmp"):
            os.remove(zfinal + ".tmp")


if __name__ == "__main__":
    main()

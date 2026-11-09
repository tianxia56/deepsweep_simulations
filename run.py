#!/usr/bin/env python3
"""Run the whole simulation pipeline from one JSON config, start to outputs.

    python3 run.py config/test_100.json            (inside the activated env)
    ./run.sh config/test_100.json                  (sets the env up first)

Stages, in order; each skips work whose completion marker already exists, so
an interrupted run is resumed by running the same command again:

  1 plan       config -> plan/*.tsv (every reference set, neutral batch and
               sweep batch, with its seeds); seed ranges checked disjoint
  2 reference  pipeline/01_make_neutral_bins.py per (model, region, target):
               neutral cosi2 sims -> per-sim pools -> exact-k bin tables
  3 neutral    pipeline/run_neutral_batch.py per batch: neutral cosi2 sims
               scored against stage 2's bins -> zipped per-sim tables
  4 sweeps     pipeline/simulate_sweep_batch.py per batch: selected cosi2
               sims -> TPED zip + par_inputs + drop list
  5 score      pipeline/02_recompute_batch.py per sweep batch against stage
               2's bins -> zipped per-sim tables; check_output.py
               (independent re-check); make_metafile.py (one causal row/sim)

Outputs go under the config's output_dir; OUTPUTS.txt lists them.
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
from datetime import datetime

import pandas as pd

REPO = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(REPO, "pipeline")
MODELS = os.path.join(REPO, "models")
SEED_MAX = 2**32                     # cosi2 truncates seeds to 32 bits
BINS_RETRY_STRIDE = 100_000_000      # 01_make_neutral_bins.SEED_RETRY_STRIDE
BINS_RETRIES_CHECKED = 9             # reseeds per replicate kept below the neutral range
NEUTRAL_SEED_STRIDE = 16             # run_neutral_batch: seed0 + sim*16 + attempt
NEUTRAL_BLOCK_MIN = 20_000           # run_neutral_batch asserts sim*16+15 < 20,000
SWEEP_SEED_STRIDE = 16               # simulate_sweep_batch: seed0 + sim*16 + attempt


def log(m):
    print(f"[{datetime.now():%F %T}] {m}", flush=True)


def die(m):
    sys.exit(f"ERROR: {m}")


# ---------------------------------------------------------------------------
# 1. plan
# ---------------------------------------------------------------------------
def par_for(model, mb):
    p = os.path.join(MODELS, f"{model}.par" if mb == 3 else f"{model}_{mb}mb.par")
    if not os.path.exists(p):
        die(f"no model file {p}")
    return p


def plan_reference(cfg):
    rb = cfg["reference_bins"]
    rows, c = [], 0
    for m in cfg["models"]:
        for mb in m["region_mb"]:
            for pop in m["target_pops"]:
                rows.append(dict(combo=f"{m['name']}_{mb}mb_pop{pop}", model=m["name"],
                                 region_mb=mb, target_pop=pop,
                                 pop_ids=",".join(map(str, m["pop_ids"])),
                                 par=par_for(m["name"], mb), n_sims=rb["n_sims"],
                                 seed0=rb["seed_base"] + c * rb["seed_step"]))
                c += 1
    if rb["n_sims"] > rb["seed_step"]:
        die("reference_bins.n_sims > seed_step: reference sets would share seeds")
    return rows


def plan_neutral(cfg, ref):
    """Same rule as the 2026oct neutral corpus (neutral_2026oct/plan_batches.py):
    per (model, Mb), n sims split as evenly as possible over the target pops,
    batches of at most max_batch, batch id {pop}{Mb}0{pop}{serial:03d}."""
    nc = cfg["neutral"]
    if nc["max_batch"] * NEUTRAL_SEED_STRIDE > NEUTRAL_BLOCK_MIN or nc["seed_block"] < NEUTRAL_BLOCK_MIN:
        die(f"neutral: need max_batch*16 <= {NEUTRAL_BLOCK_MIN} <= seed_block")
    groups = {}
    for r in ref:
        groups.setdefault((r["model"], r["region_mb"]), []).append(r)
    rows, g = [], 0
    for (model, mb), tgts in groups.items():
        n_total, k = nc["n_sims_per_model_per_mb"], len(tgts)
        per = [n_total // k + (1 if i < n_total % k else 0) for i in range(k)]
        for r, n in zip(sorted(tgts, key=lambda r: r["target_pop"]), per):
            if n == 0:
                continue
            nb = -(-n // nc["max_batch"])
            sizes = [n // nb + (1 if i < n % nb else 0) for i in range(nb)]
            par_base = os.path.basename(r["par"])[:-4]
            for s, size in enumerate(sizes):
                pop = r["target_pop"]
                batch = f"{pop}{mb}0{pop}{s:03d}"
                rows.append(dict(idx=g, model=model, region_mb=mb, target_pop=pop,
                                 pop_ids=r["pop_ids"], combo=r["combo"], par_base=par_base,
                                 batch=batch, n_sims=size,
                                 seed0=nc["seed_base"] + g * nc["seed_block"],
                                 bins_dir=os.path.join(cfg["_out"], "reference_bins", r["combo"]),
                                 zip_name=f"{par_base}.par_batch{batch}_cms_stats_all_{cfg['date']}.zip",
                                 seeds_name=f"{par_base}.par_batch{batch}_neutral_seeds_{cfg['date']}.tsv"))
                g += 1
    return rows


def plan_sweeps(cfg, ref):
    sw = cfg["sweeps"]
    combos = {(r["model"], r["region_mb"], r["target_pop"]): r for r in ref}
    rows = []
    if "batch_table" in sw:
        t = pd.read_csv(os.path.join(REPO, sw["batch_table"]), sep="\t", dtype={"batch_id": str})
        for _, b in t.iterrows():
            key = (b.model, int(b.region_mb), int(b.target_pop))
            if key not in combos:
                die(f"sweep batch {b.model}:{b.batch_id} has no reference set in this config")
            rows.append(dict(model=b.model, batch_id=b.batch_id, region_mb=int(b.region_mb),
                             target_pop=int(b.target_pop), pop_ids=b.pop_ids,
                             par=par_for(b.model, int(b.region_mb)),
                             n_sims=int(b.n_sims_planned), seed0=int(b.seed0), sweep=b.sweep))
    for d_i, d in enumerate(sw.get("designs", [])):
        names = [m["name"] for m in cfg["models"]] if d["models"] == "all" else d["models"]
        for (model, mb, pop), r in combos.items():
            if model not in names or (d.get("region_mb", "all") != "all" and mb not in d["region_mb"]) \
                    or (d.get("target_pops", "all") != "all" and pop not in d["target_pops"]):
                continue
            n = d["n_sims"]
            nb = -(-n // sw["max_batch"])
            sizes = [n // nb + (1 if i < n % nb else 0) for i in range(nb)]
            for b, size in enumerate(sizes):
                rows.append(dict(model=model, batch_id=f"{pop}{mb}{d_i:03d}{b:05d}",
                                 region_mb=mb, target_pop=pop, pop_ids=r["pop_ids"],
                                 par=r["par"], n_sims=size, seed0=None,
                                 sweep=sw["template"].format(pop=pop, age=d["age"],
                                                             daf_window=d["daf_window"])))
    nxt = sw.get("seed_base")
    for r in rows:
        if r["seed0"] is None:
            if r["n_sims"] * SWEEP_SEED_STRIDE > sw["seed_block"]:
                die("sweeps: max_batch*16 > seed_block")
            r["seed0"], nxt = nxt, nxt + sw["seed_block"]
    keys = [(r["model"], r["batch_id"]) for r in rows]
    if len(set(keys)) != len(keys):
        die("sweep batch ids repeat within a model")
    return rows


def check_seeds(ref, neut, sweeps):
    """Every seed any stage can draw lies in its own interval, below 2**32."""
    iv = []
    for r in ref:
        lo, hi = r["seed0"], r["seed0"] + r["n_sims"]
        iv += [(lo + k * BINS_RETRY_STRIDE, hi + k * BINS_RETRY_STRIDE, f"ref {r['combo']} retry {k}")
               for k in range(BINS_RETRIES_CHECKED + 1)]
    iv += [(r["seed0"], r["seed0"] + r["n_sims"] * NEUTRAL_SEED_STRIDE, f"neutral {r['model']}:{r['batch']}")
           for r in neut]
    iv += [(r["seed0"], r["seed0"] + (r["n_sims"] + 1) * SWEEP_SEED_STRIDE, f"sweep {r['model']}:{r['batch_id']}")
           for r in sweeps]
    iv.sort()
    for (a0, a1, an), (b0, b1, bn) in zip(iv, iv[1:]):
        if b0 < a1:
            die(f"seed ranges overlap: {an} [{a0}, {a1}) and {bn} [{b0}, {b1})")
    if iv and (iv[0][0] < 1 or iv[-1][1] > SEED_MAX):
        die(f"seeds outside [1, 2**32): {iv[0][2]} / {iv[-1][2]}")


def write_tsv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def plan(cfg):
    ref = plan_reference(cfg)
    neut = plan_neutral(cfg, ref)
    sweeps = plan_sweeps(cfg, ref)
    check_seeds(ref, neut, sweeps)
    d = os.path.join(cfg["_out"], "plan")
    os.makedirs(d, exist_ok=True)
    write_tsv(ref, os.path.join(d, "reference_sets.tsv"))
    if neut:
        write_tsv(neut, os.path.join(d, "neutral_batches.tsv"))
    if sweeps:
        write_tsv(sweeps, os.path.join(d, "sweep_batches.tsv"))
    log(f"plan: {len(ref)} reference sets x {cfg['reference_bins']['n_sims']} sims; "
        f"{len(neut)} neutral batches, {sum(r['n_sims'] for r in neut):,} sims; "
        f"{len(sweeps)} sweep batches, {sum(r['n_sims'] for r in sweeps):,} sims -> {d}")
    return ref, neut, sweeps


# ---------------------------------------------------------------------------
# stages 2-5
# ---------------------------------------------------------------------------
def run(cmd, logf, env=None):
    """Run one pipeline step; its output goes to logf. Non-zero exit stops the run."""
    os.makedirs(os.path.dirname(logf), exist_ok=True)
    t0 = time.time()
    with open(logf, "a") as fh:
        fh.write(f"\n### {datetime.now():%F %T}  {' '.join(cmd)}\n")
        fh.flush()
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    if rc != 0:
        tail = open(logf).read().splitlines()[-15:]
        die(f"step failed (exit {rc}) after {time.time() - t0:.0f} s: {' '.join(cmd[:3])} ...\n"
            f"  log: {logf}\n    " + "\n    ".join(tail))
    return time.time() - t0


def stage_reference(cfg, ref, env):
    for r in ref:
        out = os.path.join(cfg["_out"], "reference_bins", r["combo"])
        man = os.path.join(out, "MANIFEST.txt")
        if os.path.exists(man) and f"n_sims_asked   : {r['n_sims']}\n" in open(man).read():
            log(f"reference {r['combo']}: done")
            continue
        log(f"reference {r['combo']}: {r['n_sims']} neutral sims, seed0 {r['seed0']}")
        secs = run(["python3", os.path.join(PIPE, "01_make_neutral_bins.py"),
                    "--model", r["model"], "--par-file", r["par"],
                    "--target-pop", str(r["target_pop"]), "--pop-ids", r["pop_ids"],
                    "--region-mb", str(r["region_mb"]), "--n-sims", str(r["n_sims"]),
                    "--seed0", str(r["seed0"]), "--resume",
                    "--bins-root", os.path.join(cfg["_out"], "reference_bins"),
                    "--pools-root", os.path.join(cfg["_out"], "reference_pools", r["combo"]),
                    "--work-root", os.path.join(cfg["_tmp"], f"ref_{r['combo']}"),
                    "--workers", str(cfg["workers"])],
                   os.path.join(cfg["_out"], "logs", f"reference_{r['combo']}.log"), env)
        log(f"  -> {out} ({secs:.0f} s)")


def stage_neutral(cfg, neut, env):
    plan_f = os.path.join(cfg["_out"], "plan", "neutral_batches.tsv")
    for r in neut:
        state = os.path.join(cfg["_out"], "neutral", "state", f"{r['model']}__{r['batch']}.json")
        if os.path.exists(state):
            log(f"neutral {r['model']}:{r['batch']}: done")
            continue
        log(f"neutral {r['model']}:{r['batch']}: {r['n_sims']} sims ({r['combo']})")
        secs = run(["python3", os.path.join(PIPE, "run_neutral_batch.py"),
                    "--batches", plan_f, "--row", str(r["idx"]),
                    "--out-root", os.path.join(cfg["_out"], "neutral"),
                    "--state-root", os.path.join(cfg["_out"], "neutral", "state"),
                    "--workers", str(cfg["workers"])],
                   os.path.join(cfg["_out"], "logs", f"neutral_{r['model']}_{r['batch']}.log"), env)
        log(f"  -> {os.path.join(cfg['_out'], 'neutral', r['model'], r['zip_name'])} ({secs:.0f} s)")


def stage_sweeps(cfg, sweeps, env):
    replay = None
    if cfg["sweeps"].get("replay_sims"):
        replay = pd.read_csv(os.path.join(REPO, cfg["sweeps"]["replay_sims"]), sep="\t",
                             dtype={"batch_id": str})
    for r in sweeps:
        m, b = r["model"], r["batch_id"]
        tped_dir = os.path.join(cfg["_out"], "sweeps", "tped", m, b)
        stats_dir = os.path.join(cfg["_out"], "sweeps", "stats", m)
        meta_dir = os.path.join(cfg["_out"], "sweeps", "metafiles", m)
        done = os.path.join(cfg["_out"], "sweeps", "state", f"{m}__{b}.json")
        if os.path.exists(done):
            log(f"sweep {m}:{b}: done")
            continue
        lg = os.path.join(cfg["_out"], "logs", f"sweep_{m}_{b}.log")
        par_base = os.path.basename(r["par"])[:-4]
        cmd = ["python3", os.path.join(PIPE, "simulate_sweep_batch.py"),
               "--model", m, "--par", r["par"], "--region-mb", str(r["region_mb"]),
               "--target-pop", str(r["target_pop"]), "--pop-ids", r["pop_ids"],
               "--batch-id", b, "--n-sims", str(r["n_sims"]), "--seed0", str(r["seed0"]),
               "--sweep", r["sweep"], "--out-dir", tped_dir, "--workers", str(cfg["workers"])]
        if replay is not None:
            mine = replay[(replay.model == m) & (replay.batch_id == b)]
            rf = os.path.join(tped_dir + "_replay.tsv")
            os.makedirs(os.path.dirname(rf), exist_ok=True)
            mine[["sim_id", "seed"]].to_csv(rf, sep="\t", index=False)
            cmd += ["--replay", rf]
        log(f"sweep {m}:{b}: simulate {r['n_sims'] if replay is None else len(mine)} sims")
        run(cmd, lg, env)
        summ = json.load(open(os.path.join(tped_dir, f"{par_base}_batch{b}_summary.json")))
        if summ["to_score"] == 0:
            die(f"sweep {m}:{b}: no sim passed the causal DAF gate; nothing to score")
        log(f"  simulated {summ['simulated']}, skipped {len(summ['skipped'])}, "
            f"outside causal DAF [0.05, 0.95] {summ['dropped_for_scoring']}; scoring {summ['to_score']}")
        drop = summ["drop_list"]
        run(["python3", os.path.join(PIPE, "02_recompute_batch.py"), "--model", m,
             "--batch-id", b, "--input-dir", tped_dir, "--drop-list", drop,
             "--bins-root", os.path.join(cfg["_out"], "reference_bins"),
             "--work-root", os.path.join(cfg["_tmp"], f"score_{m}_{b}"),
             "--out-dir", stats_dir, "--workers", str(cfg["workers"])], lg, env)
        z = os.path.join(stats_dir, f"{m}_batch{b}_cms_stats_recomputed.zip")
        run(["python3", os.path.join(PIPE, "check_output.py"), m, b, tped_dir, z, drop], lg, env)
        os.makedirs(meta_dir, exist_ok=True)
        meta = os.path.join(meta_dir, f"{par_base}.par_batch{b}_par_inputs.tsv")
        run(["python3", os.path.join(PIPE, "make_metafile.py"), m, b, z,
             summ["par_inputs"], drop, meta], lg, env)
        os.makedirs(os.path.dirname(done), exist_ok=True)
        with open(done + ".tmp", "w") as fh:
            json.dump(dict(model=m, batch=b, stats_zip=z, metafile=meta,
                           tped_zip=summ["zip"], scored=summ["to_score"],
                           finished=datetime.now().isoformat(timespec="seconds")), fh, indent=1)
        os.replace(done + ".tmp", done)
        log(f"  -> {z}")
        if not cfg.get("keep_sweep_tpeds", True):
            os.remove(summ["zip"])


def write_outputs(cfg):
    out = cfg["_out"]
    lines = [f"# outputs of {cfg['run_name']}  ({datetime.now():%F %T})", ""]
    for title, pat in [("reference bin tables", "reference_bins/*/*_bin.csv"),
                       ("reference manifests", "reference_bins/*/MANIFEST.txt"),
                       ("neutral stats zips", "neutral/*/*_cms_stats_all_*.zip"),
                       ("neutral seed tables", "neutral/*/*_neutral_seeds_*.tsv"),
                       ("sweep TPED zips", "sweeps/tped/*/*/*_raw_tped.zip"),
                       ("sweep par_inputs", "sweeps/tped/*/*/*_par_inputs.tsv"),
                       ("sweep stats zips", "sweeps/stats/*/*_cms_stats_recomputed.zip"),
                       ("sweep causal metafiles", "sweeps/metafiles/*/*_par_inputs.tsv")]:
        fs = sorted(glob.glob(os.path.join(out, pat)))
        lines.append(f"## {title}: {len(fs)}")
        lines += fs + [""]
    p = os.path.join(out, "OUTPUTS.txt")
    open(p, "w").write("\n".join(lines))
    log(f"output list: {p}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--only", choices=["plan", "reference", "neutral", "sweeps"],
                    help="run one stage (its inputs must exist)")
    ap.add_argument("--workers", type=int, help="override the config's workers")
    a = ap.parse_args()

    cfg = json.load(open(a.config))
    cfg["workers"] = a.workers or cfg.get("workers") or os.cpu_count()
    od = cfg["output_dir"]
    cfg["_out"] = os.path.abspath(od if os.path.isabs(od) else os.path.join(REPO, od))
    cfg["_tmp"] = os.path.join(cfg["_out"], "tmp")
    if cfg["min_sites_per_bin"] != 5000 and not cfg.get("test"):
        die("min_sites_per_bin != 5000 is allowed only with \"test\": true")
    if cfg["reference_bins"]["n_sims"] < 1000 and cfg["min_sites_per_bin"] >= 5000:
        die(f"reference_bins.n_sims = {cfg['reference_bins']['n_sims']} cannot fill "
            f"{cfg['min_sites_per_bin']} sites per bin, so every z-score would be NaN and "
            f"scoring would stop; for a test also set \"min_sites_per_bin\": 50 and \"test\": true")
    os.makedirs(cfg["_tmp"], exist_ok=True)
    shutil.copy(a.config, os.path.join(cfg["_out"], "config.json"))
    env = dict(os.environ, TMPDIR=cfg["_tmp"], PYTHONPATH=PIPE,
               DEEPSWEEP_MIN_SITES_PER_BIN=str(cfg["min_sites_per_bin"]))
    log(f"run {cfg['run_name']} -> {cfg['_out']} ({cfg['workers']} workers)")

    ref, neut, sweeps = plan(cfg)
    if a.only == "plan":
        return
    if a.only in (None, "reference"):
        stage_reference(cfg, ref, env)
    if a.only in (None, "neutral"):
        stage_neutral(cfg, neut, env)
    if a.only in (None, "sweeps"):
        stage_sweeps(cfg, sweeps, env)
    shutil.rmtree(cfg["_tmp"], ignore_errors=True)
    write_outputs(cfg)
    log("finished")


if __name__ == "__main__":
    main()

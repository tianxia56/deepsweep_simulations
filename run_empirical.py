#!/usr/bin/env python3
"""Run the empirical pipeline from one JSON config: 1000 Genomes phase 3
VCFs (GRCh37) -> normalized per-site component statistics per population.

    python3 run_empirical.py config/empirical.json        (inside the activated env)
    ./run.sh --empirical config/empirical.json            (sets the env up first)

Stages; each skips (population, chromosome) tasks already done, so rerunning
the same command resumes:

  1 maps       convert_pyrho_maps.py    pyrho hg19 maps -> maps/<POP>/<POP>_map_chr_<N>.txt
  2 gaps       build_gap_mask.py        intervals with no map data in any core
                                        population (YRI, CEU, CHB, BEB) -> gaps.tsv
  3 ancestral  make_unified_aa_map.py   ARG ancestral allele, else the global
                                        1000 Genomes major allele -> aa/chr<N>_AA.tsv.gz
  4 tped       vcf_to_tped.py           VCF -> biallelic allele TPED per population
  5 annotate   annotate_tped.py         polarise (derived = 1), interpolate cM, drop
                                        unmapped/gap sites -> cm_tped/, cm_hmap/
  6 stats      emp_stats.py             raw statistics per (target, chromosome)
  7 bins       build_bins.py            genome-wide exact-k tables per target
  8 normalize  normalize.py             z-scores -> final/<POP>.<N>.tsv
  9 tables     stage_tables.py          one table per population
               maf_filter.py            + the maf > 0.01 version (if "maf001": true)

Stages 4-5 run for every target and every outgroup it needs.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

REPO = os.path.dirname(os.path.abspath(__file__))
EMP = os.path.join(REPO, "empirical")
sys.path.insert(0, EMP)
from emp_stats import CORE4, SUPERPOP, outgroups_for  # noqa: E402

VCF_NAME = "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"


def log(m):
    print(f"[{datetime.now():%F %T}] {m}", flush=True)


def die(m):
    sys.exit(f"ERROR: {m}")


def run(cmd, logf, env):
    os.makedirs(os.path.dirname(logf), exist_ok=True)
    with open(logf, "a") as fh:
        fh.write(f"\n### {datetime.now():%F %T}  {' '.join(cmd)}\n")
        fh.flush()
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    if rc != 0:
        tail = open(logf).read().splitlines()[-12:]
        raise RuntimeError(f"exit {rc}: {' '.join(cmd[:3])} ...\n  log: {logf}\n    "
                           + "\n    ".join(tail))


def parallel(name, tasks, workers):
    """tasks: [(label, done_marker or None, fn)]; skips tasks whose marker exists."""
    todo = [t for t in tasks if not (t[1] and os.path.exists(t[1]))]
    log(f"{name}: {len(tasks) - len(todo)} done, {len(todo)} to run")
    errors = []
    with ThreadPoolExecutor(max(1, workers)) as ex:
        futs = {ex.submit(fn): (label, marker) for label, marker, fn in todo}
        for i, f in enumerate(as_completed(futs), 1):
            label, marker = futs[f]
            try:
                f.result()
                if marker:
                    os.makedirs(os.path.dirname(marker), exist_ok=True)
                    open(marker, "w").write(datetime.now().isoformat(timespec="seconds") + "\n")
                log(f"  {name} {label} ok ({i}/{len(todo)})")
            except Exception as e:                                   # noqa: BLE001
                errors.append(f"{label}: {e}")
                log(f"  {name} {label} FAILED")
    if errors:
        die(f"{len(errors)} {name} task(s) failed:\n" + "\n".join(errors[:5]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--workers", type=int)
    a = ap.parse_args()
    cfg = json.load(open(a.config))
    W = a.workers or cfg.get("workers") or os.cpu_count()
    od = cfg["output_dir"]
    out = os.path.abspath(od if os.path.isabs(od) else os.path.join(REPO, od))
    inp = {k: os.path.abspath(v if os.path.isabs(os.path.expanduser(v)) else os.path.join(REPO, v))
           for k, v in ((k, os.path.expanduser(v)) for k, v in cfg["inputs"].items())}
    targets, chroms = cfg["targets"], [int(c) for c in cfg["chromosomes"]]
    for t in targets:
        if t not in SUPERPOP:
            die(f"target {t} is not one of the 22 scored populations ({', '.join(sorted(SUPERPOP))})")
    pops = sorted(set(targets) | {o for t in targets for o in outgroups_for(t)})
    map_pops = sorted(set(pops) | set(CORE4))
    for k in ("vcf_dir", "panel", "exclude_rsids", "ancestral_dir", "pyrho_hg19_dir"):
        if not os.path.exists(inp[k]):
            die(f"inputs.{k} not found: {inp[k]} (see empirical/download_inputs.sh)")
    for c in chroms:
        if not os.path.exists(os.path.join(inp["vcf_dir"], VCF_NAME.format(c=c))):
            die(f"missing VCF for chr{c} in {inp['vcf_dir']}")
    os.makedirs(out, exist_ok=True)
    json.dump(cfg, open(os.path.join(out, "config.json"), "w"), indent=2)
    tmp = os.path.join(out, "tmp")
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, TMPDIR=tmp, PYTHONPATH=os.path.join(REPO, "pipeline"))
    L = lambda n: os.path.join(out, "logs", n + ".log")                     # noqa: E731
    S = lambda n: os.path.join(out, "state", n + ".done")                   # noqa: E731
    py = lambda s: ["python3", os.path.join(EMP, s)]                        # noqa: E731
    log(f"empirical run {cfg['run_name']} -> {out} | targets {targets} | processed pops "
        f"{pops} | chromosomes {chroms[0]}..{chroms[-1]} ({len(chroms)}) | {W} workers")
    if len(chroms) < 22:
        log("NOTE: fewer than 22 autosomes; the genome-wide bins are built from these only (test)")

    maps, gaps, aa = (os.path.join(out, d) for d in ("maps", "gaps.tsv", "aa"))
    parallel("maps", [(p, S(f"maps_{p}"), lambda p=p: run(
        py("convert_pyrho_maps.py") + ["--pyrho-dir", inp["pyrho_hg19_dir"], "--out-dir", maps,
                                       "--pops", p], L(f"maps_{p}"), env)) for p in map_pops], W)
    parallel("gaps", [("core4", S("gaps"), lambda: run(
        py("build_gap_mask.py") + ["--pops", ",".join(CORE4), "--map-dir", maps, "--out", gaps],
        L("gaps"), env))], 1)
    parallel("ancestral", [(f"chr{c}", S(f"aa_{c}"), lambda c=c: run(
        py("make_unified_aa_map.py") + ["--chr", str(c), "--raw_vcf_dir", inp["vcf_dir"],
                                        "--ering_dir", inp["ancestral_dir"], "--out_dir", aa],
        L(f"aa_{c}"), env)) for c in chroms], W)
    parallel("tped", [(f"{p}.{c}", S(f"tped_{p}.{c}"), lambda p=p, c=c: run(
        py("vcf_to_tped.py") + ["--vcf", os.path.join(inp["vcf_dir"], VCF_NAME.format(c=c)),
                                "--panel", inp["panel"], "--pop", p,
                                "--exclude-rsids", inp["exclude_rsids"],
                                "--out-tped", os.path.join(out, "tped", f"{p}.{c}.tped")],
        L(f"tped_{p}.{c}"), env)) for p in pops for c in chroms], W)
    parallel("annotate", [(f"{p}.{c}", S(f"annotate_{p}.{c}"), lambda p=p, c=c: run(
        py("annotate_tped.py") + ["--pop", p, "--chroms", str(c),
                                  "--tped_dir", os.path.join(out, "tped"), "--rmap_dir", maps,
                                  "--gap_file", gaps, "--unified_aa_dir", aa,
                                  "--out_aa_tped_dir", os.path.join(out, "cm_tped"),
                                  "--out_hmap_dir", os.path.join(out, "cm_hmap")],
        L(f"annotate_{p}.{c}"), env)) for p in pops for c in chroms], W)
    if not cfg.get("keep_allele_tpeds", False):
        for p in pops:
            for c in chroms:
                f = os.path.join(out, "tped", f"{p}.{c}.tped")
                if os.path.exists(f) and os.path.exists(S(f"annotate_{p}.{c}")):
                    os.remove(f)
    parallel("stats", [(f"{t}.{c}", None, lambda t=t, c=c: run(
        py("emp_stats.py") + ["--target", t, "--chrom", str(c),
                              "--tped-dir", os.path.join(out, "cm_tped"),
                              "--hmap-dir", os.path.join(out, "cm_hmap"),
                              "--out-dir", os.path.join(out, "raw"),
                              "--state-dir", os.path.join(out, "state_stats")],
        L(f"stats_{t}.{c}"), env)) for t in targets for c in chroms], W)
    parallel("bins", [(t, S(f"bins_{t}"), lambda t=t: run(
        py("build_bins.py") + ["--target", t, "--raw-dir", os.path.join(out, "raw"),
                               "--out-dir", os.path.join(out, "bins", t),
                               "--chroms", ",".join(map(str, chroms))],
        L(f"bins_{t}"), env)) for t in targets], W)
    parallel("normalize", [(f"{t}.{c}", None, lambda t=t, c=c: run(
        py("normalize.py") + ["--target", t, "--chrom", str(c),
                              "--raw-dir", os.path.join(out, "raw"),
                              "--bins-dir", os.path.join(out, "bins", t),
                              "--out-dir", os.path.join(out, "final"),
                              "--state-dir", os.path.join(out, "state_norm")],
        L(f"normalize_{t}.{c}"), env)) for t in targets for c in chroms], W)
    if len(chroms) == 22:
        tab = os.path.join(out, "tables")
        parallel("tables", [(t, S(f"tables_{t}"), lambda t=t: run(
            py("stage_tables.py") + ["--pop", t, "--final-dir", os.path.join(out, "final"),
                                     "--out-dir", tab], L(f"tables_{t}"), env)) for t in targets], W)
        if cfg.get("maf001", True):
            parallel("maf001", [(t, S(f"maf001_{t}"), lambda t=t: run(
                py("maf_filter.py") + [os.path.join(tab, f"{t}.empirical.2026oct.tsv.gz"),
                                       os.path.join(out, "tables_maf001",
                                                    f"{t}.empirical.2026oct.maf001.tsv.gz")],
                L(f"maf001_{t}"), env)) for t in targets], W)
    else:
        log("tables: skipped (genome-wide tables need all 22 autosomes); per-chromosome "
            f"tables are in {os.path.join(out, 'final')}")
    log(f"finished -> {out}")


if __name__ == "__main__":
    main()

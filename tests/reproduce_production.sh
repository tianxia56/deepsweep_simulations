#!/bin/bash
# Regenerate production data with this repo and compare with what production
# wrote. Exit 0 only if every comparison passes. Needs the env (run.sh or
# env/setup_env.sh first); ~2 min on 4 cores.
#
#   bash tests/reproduce_production.sh
#
# A. Selected sweep, gravel 3 Mb target pop 2, production refill batch
#    230210020445, sim 127 (seed 3821012977):
#      simulate_sweep_batch.py --replay  -> the 3 TPEDs equal the archived ones
#                                           (md5); deri_gen, sel_gen, s, daf equal
#                                           the production metafile
#      02_recompute_batch.py against the production bins that scored it
#                                        -> every site of every column equals the
#                                           production per-sim table
# B. Reference pool, gravel_3mb_pop1, replicates 0-2 (seed0 31000000):
#      01_make_neutral_bins.py           -> per-sim (daf, value) arrays of all six
#                                           statistics equal the production pool
#                                           (XP-EHH: the keep-target-monomorphic
#                                           pool the production bins came from)
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$REPO/tests/reference"
source "$REPO/env/activate.sh"
export PYTHONPATH="$REPO/pipeline"
W=$(mktemp -d "${TMPDIR:-/tmp}/ds_repro.XXXXXX"); trap 'rm -rf "$W"' EXIT
export TMPDIR="$W"
N="${WORKERS:-4}"

echo "== A. sweep sim gravel 230210020445 #127"
S="$REF/sweep_gravel_230210020445_sim127"
printf 'sim_id\tseed\n127\t3821012977\n' > "$W/replay.tsv"
python3 "$REPO/pipeline/simulate_sweep_batch.py" --model gravel --par "$REPO/models/gravel.par" \
  --region-mb 3 --target-pop 2 --pop-ids 1,2,3 --batch-id 230210020445 --n-sims 1 \
  --seed0 3821010944 --replay "$W/replay.tsv" --out-dir "$W/tped" --workers 1 \
  --sweep 'pop_event sweep_mult_standing "sweep" 2 U(0, 2100) E(20) .5 .5-.7 2 U(0, 2100)' > "$W/sim.log"
mkdir -p "$W/bins_root/gravel_3mb_pop2"
( cd "$W" && unzip -q "$S/production_stats.zip" 'bins/*' && mv bins/* bins_root/gravel_3mb_pop2/ )
python3 "$REPO/pipeline/02_recompute_batch.py" --model gravel --batch-id 230210020445 \
  --input-dir "$W/tped" --drop-list "$W/tped/gravel_batch230210020445_drop.tsv" \
  --bins-root "$W/bins_root" --work-root "$W/work" --out-dir "$W/stats" --workers "$N" > "$W/score.log" 2>&1 \
  || { tail -20 "$W/score.log"; exit 1; }

echo "== B. reference pool gravel_3mb_pop1 replicates 0-2"
python3 "$REPO/pipeline/01_make_neutral_bins.py" --model gravel --par-file "$REPO/models/gravel.par" \
  --target-pop 1 --pop-ids 1,2,3 --region-mb 3 --n-sims 3 --seed0 31000000 \
  --bins-root "$W/bins_new" --pools-root "$W/pools_new" --work-root "$W/simwork" \
  --workers "$N" > "$W/pool.log" 2>&1 || { tail -20 "$W/pool.log"; exit 1; }

python3 - "$W" "$REF" "$REPO" <<'PY'
import glob, hashlib, io, os, sys, zipfile
import numpy as np, pandas as pd
w, ref, repo = sys.argv[1:4]
fail = []

# A1: TPEDs byte-identical to the archived ones
want = dict(l.split()[::-1] for l in open(f"{repo}/env/reference_md5.txt") if l.strip() and not l.startswith("#"))
with zipfile.ZipFile(f"{w}/tped/gravel_batch230210020445_raw_tped.zip") as z:
    for p in (1, 2, 3):
        got = hashlib.md5(z.read(f"sel/sel.hap.127_0_{p}.tped")).hexdigest()
        ok = got == want[f"sweep_gravel_3mb_pop2_seed3821012977_pop{p}.tped"]
        print(f"A1 TPED pop {p}: {'identical' if ok else 'DIFFERENT'}"); ok or fail.append(f"TPED pop {p}")
# A2: sampled sweep parameters and causal DAF equal the production metafile
pi = pd.read_csv(glob.glob(f"{w}/tped/*_par_inputs.tsv")[0], sep="\t").iloc[0]
pm = pd.read_csv(f"{ref}/sweep_gravel_230210020445_sim127/production_metafile.tsv", sep="\t").iloc[0]
for c in ("daf", "deri_gen", "sel_gen", "s"):
    ok = np.isclose(float(pi[c]), float(pm[c]), atol=1e-6, rtol=0)
    print(f"A2 {c}: {pi[c]} vs production {pm[c]}: {'equal' if ok else 'DIFFERENT'}"); ok or fail.append(c)
# A3: every site of every column equals the production per-sim table
name = "gravel_batch230210020445_cms_stats_127.tsv"
new = pd.read_csv(io.BytesIO(zipfile.ZipFile(f"{w}/stats/gravel_batch230210020445_cms_stats_recomputed.zip").read(name)), sep="\t")
old = pd.read_csv(io.BytesIO(zipfile.ZipFile(f"{ref}/sweep_gravel_230210020445_sim127/production_stats.zip").read(name)), sep="\t")
if list(new.columns) != list(old.columns) or len(new) != len(old):
    fail.append("table shape"); print(f"A3 shape {new.shape} vs {old.shape}")
else:
    worst = 0.0
    for c in new.columns:
        a, b = new[c].to_numpy(float), old[c].to_numpy(float)
        same_na = (np.isnan(a) == np.isnan(b)).all()
        d = np.nanmax(np.abs(a - b)) if (~np.isnan(a)).any() else 0.0
        worst = max(worst, d)
        if not same_na or d > 0:
            fail.append(f"column {c}"); print(f"A3 {c}: NA pattern equal {same_na}, max |diff| {d:g}")
    print(f"A3 {len(new):,} sites x {len(new.columns)} columns: max |new - production| = {worst:g}")
# B: pooled (daf, value) arrays equal production
for f in sorted(glob.glob(f"{ref}/pools_gravel_3mb_pop1/*/*.npz")):
    stat, sim = f.split("/")[-2], os.path.basename(f)
    a, b = np.load(f), np.load(f"{w}/pools_new/{stat}/{sim}")
    if "k" in a.files:
        # XP-EHH: the production bins were rebuilt (2026-10-16, keeping
        # target-monomorphic sites) from pools keyed on the exact count k
        k_new = np.rint(b["daf"] * int(b["n_hap"])).astype(int)
        ok = np.array_equal(a["k"], k_new) and np.array_equal(a["val"], b["val"], equal_nan=True)
    else:
        ok = all(np.array_equal(a[k], b[k], equal_nan=True) for k in ("daf", "val"))
    if not ok:
        fail.append(f"pool {stat}/{sim}")
    print(f"B  pool {stat:13s} {sim}: {'identical' if ok else 'DIFFERENT'} ({len(a['val']):,} sites)")
print("REPRODUCTION OK" if not fail else f"REPRODUCTION FAILED: {fail}")
sys.exit(1 if fail else 0)
PY

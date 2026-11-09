#!/usr/bin/env python3
"""Independent check of one recomputed batch before its DONE marker is written.

    check_output.py <model> <batch_id> <input_dir> <output_zip> <drop_list>

Expected sims = sims with TPEDs in the input zip, minus the drop list's rows
for (model, batch). Shares no code with 02_recompute_batch.py on purpose.
Fails (exit 1) unless:
  - the output zip's per-sim TSV ids equal the expected set exactly
  - every TSV has exactly one selpos == 1 row, at pos == selpos, with
    causal DAF in [0.05, 0.95] (the causal MAF gate)
"""
import glob, io, re, sys, zipfile
import pandas as pd

model, batch, inp, out, drop_f = sys.argv[1:6]
selpos = int(batch[1]) * 1_000_000 // 2
zins = glob.glob(f"{inp}/*batch{batch}*tped*.zip")
assert len(zins) == 1, f"input zips: {zins}"
with zipfile.ZipFile(zins[0]) as zf:
    have = {m.group(1) for n in zf.namelist()
            if (m := re.search(r"sel\.hap\.(\d+)_0_\d+\.tped$", n))}
dl = pd.read_csv(drop_f, sep="\t", dtype={"sim_batch_no": str})
drop = set(dl[(dl.model == model) & (dl.sim_batch_no == batch)].sim_id.astype(int).astype(str))
expected = have - drop

probs = []
with zipfile.ZipFile(out) as zf:
    tsvs = {m.group(1): n for n in zf.namelist()
            if (m := re.search(rf"^{re.escape(model)}_batch{batch}_cms_stats_(\d+)\.tsv$", n))}
    got = set(tsvs)
    if got != expected:
        probs.append(f"sim set: {len(got)} written vs {len(expected)} expected; "
                     f"missing {sorted(expected - got, key=int)[:10]}, "
                     f"unexpected {sorted(got - expected, key=int)[:10]}")
    for sid, n in tsvs.items():
        d = pd.read_csv(io.BytesIO(zf.read(n)), sep="\t", usecols=["pos", "daf", "selpos"])
        c = d[d.selpos == 1]
        if len(c) != 1 or int(c.pos.iloc[0]) != selpos:
            probs.append(f"sim {sid}: {len(c)} causal rows"); continue
        if not 0.05 <= round(float(c.daf.iloc[0]), 6) <= 0.95:
            probs.append(f"sim {sid}: causal DAF {float(c.daf.iloc[0]):.3f} outside [0.05, 0.95]")

print(f"{model}:{batch}  input {len(have)}  dropped {len(have & drop)}  "
      f"expected {len(expected)}  written {len(got)}  problems {len(probs)}")
for p in probs[:20]:
    print("  ", p)
sys.exit(1 if probs else 0)

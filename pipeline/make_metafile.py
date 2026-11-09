#!/usr/bin/env python3
"""Causal-site metafile for one scored sweep batch.

One row per scored sim: the selpos == 1 row of its per-sim TSV, in the
production metafile column order. Verified against the pre-scoring
par_inputs (hard error on any mismatch):
  sim set == par_inputs sims - drop list
  daf, deri_gen, sel_gen, s equal within 1e-6

  make_metafile.py <model> <batch_id> <stats_zip> <par_inputs.tsv> <drop.tsv> <out.tsv>
"""
import json
import os
import re
import sys
import zipfile

import numpy as np
import pandas as pd

COLS = ["sim_batch_no", "sim_id", "pos", "mean_fst", "deldaf", "daf", "maf",
        "norm_nsl", "norm_ihh12", "norm_ihs_binom", "norm_delihh_binom",
        "max_xpehh", "deri_gen", "sel_gen", "s", "selpos"]
CHECK = ["daf", "deri_gen", "sel_gen", "s"]


def main():
    m, b, z, src, drop_f, out = sys.argv[1:7]
    selpos = int(b[1]) * 1_000_000 // 2
    rows = []
    with zipfile.ZipFile(z) as zf:
        for n in zf.namelist():
            if not re.fullmatch(rf"{re.escape(m)}_batch{b}_cms_stats_\d+\.tsv", n):
                continue
            with zf.open(n) as fh:
                d = pd.read_csv(fh, sep="\t", na_values="NA", dtype={"sim_batch_no": str})
            c = d[d.selpos == 1]
            assert len(c) == 1 and int(c.pos.iloc[0]) == selpos, (m, b, n, len(c))
            rows.append(c)
    meta = pd.concat(rows)[COLS].sort_values("sim_id").reset_index(drop=True)
    assert (meta.sim_batch_no == b).all()

    s = pd.read_csv(src, sep="\t")
    s = s.set_index(s.sim_id.astype(int))
    dl = pd.read_csv(drop_f, sep="\t", dtype={"sim_batch_no": str})
    drop = set(dl[(dl.model == m) & (dl.sim_batch_no == b)].sim_id.astype(int))
    want, got = set(s.index) - drop, set(meta.sim_id.astype(int))
    assert got == want, (m, b, "sim set", len(got), len(want), sorted(got ^ want)[:10])
    x = meta.set_index(meta.sim_id.astype(int))
    for c in CHECK:
        a, o = x[c].astype(float), s.loc[x.index, c].astype(float)
        bad = int((~np.isclose(a, o, atol=1e-6, rtol=0)).sum())
        assert bad == 0, (m, b, c, bad)

    for c in ["sim_id", "pos", "deri_gen", "sel_gen", "selpos"]:
        meta[c] = meta[c].astype("Int64")
    meta.to_csv(out + ".tmp", sep="\t", index=False, na_rep="NA")
    os.replace(out + ".tmp", out)
    na = {c: int(meta[c].isna().sum()) for c in COLS if meta[c].isna().any()}
    with open(out + ".check.json", "w") as fh:
        json.dump(dict(model=m, batch=b, n_sims=len(meta), dropped=len(drop),
                       source=os.path.basename(src), na_at_causal=na), fh)
    print(f"{m}:{b} {len(meta)} sims -> {out}" + (f"  NA at causal {na}" if na else ""))


if __name__ == "__main__":
    main()

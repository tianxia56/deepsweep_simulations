#!/usr/bin/env python3
"""Stage one gzipped table per population for upload, and verify it.

final/ holds 484 per-chromosome tables; the bucket layout is one file per
population, matching how the previous empirical buckets were organised.

Verification is done on the STAGED file, not on the inputs: the concatenation
is the step that can go wrong (a lost header line, a chromosome written twice,
a truncated gzip), so it is the step that is checked.

  V1 one header line, exactly the expected columns
  V2 data line count == sum of input data lines
  V3 chromosomes appear once each, in order 1..22
  V4 the gzip decompresses cleanly end to end (implied by V2/V3, which read
     the whole file)

-> <out-dir>/<POP>.empirical.2026oct.tsv.gz  + <POP>.json
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil

HEADER = ("chr\tpos\tID\tmean_fst\tdeldaf\tdaf\tmaf\tnorm_nsl\tnorm_ihh12\t"
          "norm_ihs_binom\tnorm_delihh_binom\tmax_xpehh")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", required=True)
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--level", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    out = f"{a.out_dir}/{a.pop}.empirical.2026oct.tsv.gz"
    tmp = out + ".tmp"

    expect_rows = 0
    with gzip.open(tmp, "wt", compresslevel=a.level) as w:
        w.write(HEADER + "\n")
        for c in range(1, 23):
            f = f"{a.final_dir}/{a.pop}.{c}.tsv"
            with open(f) as fh:
                head = fh.readline().rstrip("\n")
                if head != HEADER:
                    raise SystemExit(f"ERROR: {f} header mismatch")
                for line in fh:
                    w.write(line)
                    expect_rows += 1
    os.replace(tmp, out)

    # ---- verify the staged file by reading it back ----------------------
    seen, order, nrow = {}, [], 0
    with gzip.open(out, "rt") as fh:
        if fh.readline().rstrip("\n") != HEADER:
            raise SystemExit("ERROR: staged header mismatch")
        for line in fh:
            nrow += 1
            c = line[:line.index("\t")]
            if c not in seen:
                seen[c] = 1
                order.append(int(c))
            else:
                seen[c] += 1
    problems = []
    if nrow != expect_rows:
        problems.append(f"V2 rows {nrow} != {expect_rows}")
    if order != list(range(1, 23)):
        problems.append(f"V3 chromosome order {order}")
    h = hashlib.md5()
    with open(out, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    info = dict(pop=a.pop, file=os.path.basename(out),
                bytes=os.path.getsize(out), rows=nrow,
                rows_expected=expect_rows, chrom_order_ok=(order == list(range(1, 23))),
                per_chrom={str(k): v for k, v in seen.items()},
                md5=h.hexdigest(), problems=problems,
                verdict="PASS" if not problems else "FAIL")
    with open(f"{a.out_dir}/{a.pop}.json", "w") as fh:
        json.dump(info, fh, indent=1)
    print(f"{a.pop}: {info['verdict']} {nrow:,} rows, "
          f"{info['bytes']/1e9:.2f} GB {problems if problems else ''}")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())

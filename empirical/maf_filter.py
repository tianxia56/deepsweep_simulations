#!/usr/bin/env python3
"""Keep rows with maf > 0.01 (strict) of one genome-wide population table.

Python port of the production maf001/filter.sh awk: same rule, same header,
column 7 must be maf. The training corpus is NOT filtered; this filter is
applied to the empirical tables only.

  maf_filter.py <in .tsv.gz> <out .tsv.gz>
"""
import gzip
import os
import sys


def main():
    src, out = sys.argv[1:3]
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    n_in = n_out = ties = 0
    with gzip.open(src, "rt") as fh, gzip.open(out + ".part", "wt", compresslevel=6) as o:
        head = fh.readline()
        if head.split("\t")[6] != "maf":
            sys.exit("ERROR: maf is not column 7")
        o.write(head)
        for line in fh:
            n_in += 1
            m = float(line.split("\t", 7)[6])
            if m > 0.01:
                o.write(line)
                n_out += 1
            elif m == 0.01:
                ties += 1
    os.replace(out + ".part", out)
    print(f"{os.path.basename(src)}: {n_in:,} rows in, {n_out:,} kept (maf > 0.01), "
          f"{ties:,} at exactly 0.01 dropped -> {out}")


if __name__ == "__main__":
    main()

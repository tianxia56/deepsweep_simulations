#!/usr/bin/env python3
"""pyrho hg19 maps (HapMap format) -> the per-population map files the
empirical pipeline reads.

  in : <pyrho-dir>/<POP>/<POP>_recombination_map_hapmap_format_hg19_chr_<N>.txt
       Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)      (chr22 16050114 ...)
  out: <out-dir>/<POP>/<POP>_map_chr_<N>.txt
       pos  chr  cM                                        (16050114 22 ...)

Values are copied as text (no float round trip), so the output is byte-identical
to the maps used for the 2026oct corpus. Source: Spence & Song 2019, Zenodo
10.5281/zenodo.11437540 (hg19_maps.tar.gz).
"""
import argparse
import os
import sys


def convert(src, dst):
    with open(src) as fh, open(dst + ".tmp", "w") as out:
        head = fh.readline().split()
        if head != ["Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)"]:
            sys.exit(f"ERROR: {src}: unexpected header {head}")
        out.write("pos\tchr\tcM\n")
        n = 0
        for line in fh:
            c, pos, _, cm = line.rstrip("\n").split("\t")
            out.write(f"{pos}\t{c[3:] if c.startswith('chr') else c}\t{cm}\n")
            n += 1
    os.replace(dst + ".tmp", dst)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyrho-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--pops", required=True, help="comma list")
    a = ap.parse_args()
    for pop in a.pops.split(","):
        os.makedirs(os.path.join(a.out_dir, pop), exist_ok=True)
        for c in range(1, 23):
            src = os.path.join(a.pyrho_dir, pop,
                               f"{pop}_recombination_map_hapmap_format_hg19_chr_{c}.txt")
            if not os.path.exists(src):
                sys.exit(f"ERROR: map not found: {src}")
            convert(src, os.path.join(a.out_dir, pop, f"{pop}_map_chr_{c}.txt"))
        print(f"{pop}: 22 maps -> {os.path.join(a.out_dir, pop)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# File: annotate_tped.py (production name: unified_annotation_pipeline.py)
import os
import sys
import bisect
import gzip
import logging
import argparse

VALID_BASES = {'A', 'C', 'G', 'T'}
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_rmap(rmap_file):
    """-> (positions, cM floats) ascending and deduplicated.

    A missing or unusable map is a HARD ERROR. It used to return (None, None),
    after which every site on the chromosome silently received cM 0 -- a whole
    chromosome of zero genetic distance, which no downstream stage could
    detect.
    """
    if not os.path.exists(rmap_file):
        sys.exit(f"ERROR: recombination map not found: {rmap_file}")
    is_gz = rmap_file.endswith(".gz")
    open_fn = gzip.open if is_gz else open
    rows = []
    with open_fn(rmap_file, 'rt', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3 or not parts[0].isdigit():
                continue                      # header or blank
            try:
                rows.append((int(parts[0]), float(parts[2])))
            except ValueError:
                continue
    if len(rows) < 2:
        sys.exit(f"ERROR: map {rmap_file} has {len(rows)} usable rows")
    rows.sort()
    positions, cm_values = [], []
    for pos, cm in rows:
        if positions and pos == positions[-1]:
            continue                          # duplicate position: keep first
        positions.append(pos)
        cm_values.append(cm)
    if any(cm_values[i] > cm_values[i + 1] for i in range(len(cm_values) - 1)):
        sys.exit(f"ERROR: map {rmap_file} is not monotone in cM after sorting")
    return positions, cm_values


def load_gaps(path, chrom):
    """Intervals unmappable in EVERY population, dropped identically for all."""
    if not path:
        return []
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.split()
            if parts and int(parts[0]) == chrom:
                out.append((int(parts[1]), int(parts[2])))
    return sorted(out)


def in_gap(gap_starts, gaps, pos):
    i = bisect.bisect_right(gap_starts, pos) - 1 if gaps else -1
    return i >= 0 and gaps[i][0] < pos < gaps[i][1]


def cm_at(map_pos, map_cm, pos, last_i):
    """cM by LINEAR INTERPOLATION between the flanking map rows.

    The previous code took map_cm[bisect_right(map_pos, pos) - 1], the value of
    the nearest row at or below the site. That is a step function: on YRI chr22
    its 1,054,119 variants carried only 23,544 distinct cM values, so ~45
    consecutive SNPs sat at an identical genetic position. iHS, dihh, iHH12 and
    XP-EHH all integrate EHH over genetic distance, so those runs contributed
    ZERO to the integral.
    """
    i = bisect.bisect_right(map_pos, pos) - 1
    if i >= last_i:
        return map_cm[last_i]
    lo, hi = map_pos[i], map_pos[i + 1]
    return map_cm[i] + (pos - lo) * (map_cm[i + 1] - map_cm[i]) / (hi - lo)

def load_unified_aa(aa_file_gz):
    aa_map = {}
    if not os.path.exists(aa_file_gz):
        return None
    with gzip.open(aa_file_gz, 'rt', encoding='utf-8') as f:
        f.readline() # Header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[0].isdigit():
                aa_map[int(parts[0])] = parts[1]
    return aa_map

def process_chromosome(pop, chrom, in_tped, rmap_file, aa_file, out_aa_tped,
                       out_map, out_hap, gaps=()):
    logging.info(f"Pop: {pop} | Chr: {chrom}")
    if not os.path.exists(in_tped):
        logging.error(f"Input TPED missing: {in_tped}")
        return False

    map_pos, map_cm = load_rmap(rmap_file)
    aa_map = load_unified_aa(aa_file)
    if aa_map is None:
        logging.error(f"Unified AA reference missing for Chr {chrom}: {aa_file}")
        return False

    records = []
    skipped_cnt = 0
    dropped = {"before_map": 0, "after_map": 0, "in_gap": 0}
    last_i = len(map_pos) - 1
    gap_starts = [g[0] for g in gaps]
    seen_ids = set()
    renamed_dot = 0

    with open(in_tped, 'r', encoding='utf-8', buffering=1024*1024) as infile:
        for line in infile:
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) < 5:
                continue

            tped_chr, rsid, _, pos_str = fields[:4]
            genotypes = fields[4:]
            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # Assign cM from the recombination map, by interpolation. A site
            # outside the mapped range is DROPPED: there is no defensible
            # genetic position for it, and the old code gave those cM 0 (before
            # the map) or the final map value (after it), both silently.
            if pos < map_pos[0]:
                dropped["before_map"] += 1
                continue
            if pos > map_pos[last_i]:
                dropped["after_map"] += 1
                continue
            if in_gap(gap_starts, gaps, pos):
                dropped["in_gap"] += 1
                continue
            cm_val = f"{cm_at(map_pos, map_cm, pos, last_i):.10f}"

            # Ancestral allele lookup (guaranteed identical across all populations)
            aa = aa_map.get(pos)
            if not aa:
                skipped_cnt += 1
                continue

            # Convert genotypes: 0 = ancestral, 1 = derived
            valid = True
            aa_genotypes = []
            for g in genotypes:
                gu = g.upper()
                if gu == aa:
                    aa_genotypes.append('0')
                elif gu in VALID_BASES:
                    aa_genotypes.append('1')
                else:
                    valid = False
                    break

            if not valid:
                skipped_cnt += 1
                continue

            if rsid == ".":
                rsid = f"rspos{pos_str}"
                renamed_dot += 1
            if rsid in seen_ids:
                sys.exit(f"ERROR: duplicate variant ID {rsid} at {chrom}:{pos}; "
                         f"IDs must be unique because the statistics are "
                         f"joined on them")
            seen_ids.add(rsid)

            hap_str = " ".join(aa_genotypes)
            tped_line = f"{tped_chr}\t{rsid}\t{cm_val}\t{pos_str}\t{hap_str}\n"
            map_line = f"{tped_chr} {rsid} {cm_val} {pos_str}\n"
            hap_line = f"{hap_str}\n"

            records.append((pos, tped_line, map_line, hap_line))

    records.sort(key=lambda x: x[0])

    with open(out_aa_tped, 'w', encoding='utf-8', buffering=1024*1024) as f_tped, \
         open(out_map, 'w', encoding='utf-8', buffering=1024*1024) as f_map, \
         open(out_hap, 'w', encoding='utf-8', buffering=1024*1024) as f_hap:
        for _, tped_line, map_line, hap_line in records:
            f_tped.write(tped_line)
            f_map.write(map_line)
            f_hap.write(hap_line)

    logging.info(
        f"  Chr {chrom} Done: {len(records):,} written | {skipped_cnt:,} "
        f"skipped (no AA / invalid GT) | dropped before_map "
        f"{dropped['before_map']:,}, after_map {dropped['after_map']:,}, "
        f"in_gap {dropped['in_gap']:,} | renamed '.' ids {renamed_dot:,}")
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pop', required=True)
    parser.add_argument('--tped_dir', required=True)
    parser.add_argument('--rmap_dir', required=True,
                        help='selscan maps: <dir>/<pop>/<pop>_map_chr_<N>.txt')
    parser.add_argument('--chroms', default=",".join(str(c) for c in range(1, 23)))
    parser.add_argument('--gap_file', default=None,
                        help='build_gap_mask.py output: intervals unmappable '
                             'in EVERY population, dropped identically for all')
    parser.add_argument('--unified_aa_dir', required=True)
    parser.add_argument('--out_aa_tped_dir', required=True)
    parser.add_argument('--out_hmap_dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_aa_tped_dir, exist_ok=True)
    os.makedirs(args.out_hmap_dir, exist_ok=True)

    success = True
    for chrom in [int(c) for c in args.chroms.split(",")]:
        in_tped = os.path.join(args.tped_dir, f"{args.pop}.{chrom}.tped")
        rmap_file = os.path.join(args.rmap_dir, args.pop, f"{args.pop}_map_chr_{chrom}.txt")
        aa_file = os.path.join(args.unified_aa_dir, f"chr{chrom}_AA.tsv.gz")
        out_aa_tped = os.path.join(args.out_aa_tped_dir, f"{args.pop}.{chrom}.tped")
        out_map = os.path.join(args.out_hmap_dir, f"{args.pop}.{chrom}.map")
        out_hap = os.path.join(args.out_hmap_dir, f"{args.pop}.{chrom}.hap")

        gaps = load_gaps(args.gap_file, chrom)
        ok = process_chromosome(args.pop, chrom, in_tped, rmap_file, aa_file,
                                out_aa_tped, out_map, out_hap, gaps)
        if not ok:
            success = False

    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()

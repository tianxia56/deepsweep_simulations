#!/usr/bin/env python3
# File: make_unified_aa_map.py
import os
import sys
import gzip
import logging
import argparse

VALID_BASES = {'A', 'C', 'G', 'T'}
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_ering_aa(ering_file):
    aa_map = {}
    if not os.path.exists(ering_file):
        sys.exit(f"ERROR: ancestral-allele table not found: {ering_file}")
    with gzip.open(ering_file, 'rt', encoding='utf-8') as f:
        f.readline() # Header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[1].isdigit():
                al = parts[2].upper()
                if al in VALID_BASES:
                    aa_map[int(parts[1])] = al
    return aa_map

def extract_global_major_from_vcf(vcf_gz, ering_map):
    """
    Scans the raw 1KGP VCF. Only checks the first 8 columns (fast!).
    If pos is in ErinG, keeps ErinG. If missing, uses global AF to define major.
    """
    unified_map = {}
    ering_hits = 0
    fallback_hits = 0

    with gzip.open(vcf_gz, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t', 8)
            if len(parts) < 8:
                continue

            pos_str, ref, alt, info = parts[1], parts[3].upper(), parts[4].upper(), parts[7]
            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # 1. Use ErinG ARG AA if present
            if pos in ering_map:
                unified_map[pos] = ering_map[pos]
                ering_hits += 1
                continue

            # 2. Fallback to Global Major Allele from 1KGP cohort AF
            # AF is alternate allele frequency across all 2504 samples
            af = None
            for item in info.split(';'):
                if item.startswith('AF='):
                    try:
                        af = float(item.split('=')[1].split(',')[0])
                    except ValueError:
                        pass
                    break

            if af is not None and ref in VALID_BASES and alt in VALID_BASES:
                # If AF <= 0.5, REF is global major; if AF > 0.5, ALT is global major
                major_al = ref if af <= 0.5 else alt
                unified_map[pos] = major_al
                fallback_hits += 1

    logging.info(f"  ErinG hits: {ering_hits:,} | Global Major fallbacks: {fallback_hits:,}")
    return unified_map

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chr', type=int, required=True)
    parser.add_argument('--raw_vcf_dir', required=True)
    parser.add_argument('--ering_dir', required=True,
                        help="chr<N>_ErinG_AA_hg19.tsv.gz: CHR POS ErinG_ARG_AA ...")
    parser.add_argument('--out_dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    out_file = os.path.join(args.out_dir, f"chr{args.chr}_AA.tsv.gz")

    ering_file = os.path.join(args.ering_dir, f"chr{args.chr}_ErinG_AA_hg19.tsv.gz")
    raw_vcf = os.path.join(args.raw_vcf_dir, f"ALL.chr{args.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

    logging.info(f"Processing Chr {args.chr}...")
    if not os.path.exists(raw_vcf):
        sys.exit(f"ERROR: VCF not found: {raw_vcf}")
    ering_map = load_ering_aa(ering_file)
    unified_map = extract_global_major_from_vcf(raw_vcf, ering_map)

    # Write unified map
    with gzip.open(out_file, 'wt', encoding='utf-8') as out_f:
        out_f.write("POS\tAA\n")
        for pos in sorted(unified_map.keys()):
            out_f.write(f"{pos}\t{unified_map[pos]}\n")

    logging.info(f"Wrote {len(unified_map):,} variants to {out_file}")

if __name__ == '__main__':
    main()

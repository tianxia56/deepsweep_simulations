#!/usr/bin/env python3
"""One population, one chromosome: 1000 Genomes phase 3 VCF -> allele TPED.

The per-site logic is the production streamer embedded in
00run_vcf_to_tped_pipeline.sh (2026-09), unchanged, so the output reproduces
the TPEDs the 2026oct empirical corpus was built from:

  vcftools --gzvcf <vcf> --keep <samples> --min-alleles 2 --max-alleles 2
           --recode --stdout            (biallelic records, NO frequency filter:
                                         sites monomorphic in this population stay)
  0  drop records whose rsID is in the exclusion list (merged, retired rsIDs)
  1  rename IDs that do not start with "rs" to rspos<pos> ("." is left as is;
     annotate_tped.py renames it)
  2  drop a record whose ID was already seen
  3  drop a record at an already-seen position (the first record is kept)
  4  drop a record with any missing, non-0/1 or malformed genotype
  -> "<chr> <id> 0 <pos> <allele> <allele> ..." two alleles per sample, as
     REF/ALT letters; annotate_tped.py polarises them

Differences from production, none of which changes a written line: a missing
VCF, sample list or exclusion list, or a non-zero vcftools exit, is an ERROR
(production skipped or continued); vcftools stderr goes to a log file; output
is written to a temporary name and renamed only on success.

  vcf_to_tped.py --vcf ALL.chr22....vcf.gz --panel integrated_call_samples_v3.20130502.ALL.panel
                 --pop YRI --exclude-rsids exclude_rsid.txt --out-tped tped/YRI.22.tped
"""
import argparse
import os
import subprocess
import sys


def sample_list(panel, pop, out):
    """Panel rows (sample, pop, ...) of one population, in panel order."""
    with open(panel) as fh, open(out, "w") as o:
        n = 0
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0 or len(f) < 2:
                continue
            if f[1] == pop:
                o.write(f[0] + "\n")
                n += 1
    if n == 0:
        sys.exit(f"ERROR: no samples of {pop} in {panel}")
    return n


def run_stream(chrom, raw_vcf, keep_list, exclude_file, out_tped, stats_out, err_log):
    excluded_rsids = set()
    with open(exclude_file, 'r') as f:
        for line in f:
            s = line.strip()
            if s.startswith("rs"):
                excluded_rsids.add(s)

    vcftools_cmd = [
        "vcftools",
        "--gzvcf", raw_vcf,
        "--keep", keep_list,
        "--min-alleles", "2",
        "--max-alleles", "2",
        "--recode",
        "--stdout"
    ]

    with open(err_log, "w") as err:
        p_vcf = subprocess.Popen(vcftools_cmd, stdout=subprocess.PIPE, stderr=err,
                                 text=True, bufsize=1024*1024)

        initial_var = 0
        outdated_ids = 0
        ids_renamed = 0
        dup_ids = 0
        dup_pos = 0
        inv_gt = 0
        final_var = 0

        seen_ids = set()
        seen_pos = set()

        gt_idx = -1
        header_found = False

        with open(out_tped, 'w', buffering=1024*1024) as out_f:
            for line in p_vcf.stdout:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    header_found = True
                    continue
                if not header_found:
                    continue

                initial_var += 1
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue

                vcf_chr, pos, rsid_orig, ref, alt_field, _, _, _, fmt = fields[:9]
                samples = fields[9:]
                alt = alt_field.split(',')[0]

                # 0. Exclude outdated rsID
                if rsid_orig != '.' and rsid_orig in excluded_rsids:
                    outdated_ids += 1
                    continue

                # Determine GT index once
                if gt_idx == -1:
                    try:
                        gt_idx = fmt.split(':').index('GT')
                    except ValueError:
                        continue

                # 1. Rename non-rs ID
                cur_id = rsid_orig
                if cur_id != '.' and not cur_id.startswith("rs"):
                    new_id = f"rspos{pos}"
                    if new_id != cur_id:
                        cur_id = new_id
                        ids_renamed += 1

                # 2. Duplicate ID Check
                if cur_id != '.' and cur_id in seen_ids:
                    dup_ids += 1
                    continue
                if cur_id != '.':
                    seen_ids.add(cur_id)

                # 3. Duplicate Position Check
                pos_key = (vcf_chr, pos)
                if pos_key in seen_pos:
                    dup_pos += 1
                    if cur_id != '.' and cur_id in seen_ids:
                        seen_ids.remove(cur_id)
                    continue
                seen_pos.add(pos_key)

                # 4. Parse & Validate Genotypes directly to TPED alleles
                # 0 -> REF, 1 -> ALT
                valid = True
                tped_alleles = []

                for s in samples:
                    if ':' in s:
                        s_parts = s.split(':')
                        gt_str = s_parts[gt_idx] if gt_idx < len(s_parts) else '.'
                    else:
                        gt_str = s

                    if '.' in gt_str or len(gt_str) < 3:
                        valid = False
                        break

                    a1 = gt_str[0]
                    a2 = gt_str[2]

                    if a1 not in ('0', '1') or a2 not in ('0', '1'):
                        valid = False
                        break

                    tped_alleles.append(ref if a1 == '0' else alt)
                    tped_alleles.append(ref if a2 == '0' else alt)

                if not valid:
                    inv_gt += 1
                    if cur_id != '.' and cur_id in seen_ids:
                        seen_ids.remove(cur_id)
                    seen_pos.remove(pos_key)
                    continue

                # Write standard TPED format: CHR SNP_ID 0 POS A1 A2 A3 A4...
                out_f.write(f"{vcf_chr} {cur_id} 0 {pos} {' '.join(tped_alleles)}\n")
                final_var += 1

        rc = p_vcf.wait()
    if rc != 0:
        sys.exit(f"ERROR: vcftools exited {rc} on {raw_vcf}; see {err_log}")
    if not header_found:
        sys.exit(f"ERROR: vcftools wrote no #CHROM header for {raw_vcf}; see {err_log}")

    with open(stats_out, 'w') as sf:
        sf.write("chr\tinitial_var\toutdated_id\tid_renamed\tdup_id\tdup_pos\tinvalid_gt\tfinal_var\n")
        sf.write(f"{chrom}\t{initial_var}\t{outdated_ids}\t{ids_renamed}\t{dup_ids}\t{dup_pos}\t{inv_gt}\t{final_var}\n")
    return final_var


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--panel", required=True, help="integrated_call_samples_v3.20130502.ALL.panel")
    ap.add_argument("--pop", required=True)
    ap.add_argument("--exclude-rsids", required=True, help="one rsID per line (dbSNP RsMergeArch)")
    ap.add_argument("--out-tped", required=True)
    a = ap.parse_args()
    for f in (a.vcf, a.panel, a.exclude_rsids):
        if not os.path.exists(f):
            sys.exit(f"ERROR: not found: {f}")
    chrom = os.path.basename(a.vcf).split(".")[1].replace("chr", "")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_tped)), exist_ok=True)
    tmp = a.out_tped + ".tmp"
    keep = a.out_tped + ".samples.txt"
    n = sample_list(a.panel, a.pop, keep)
    try:
        rows = run_stream(chrom, a.vcf, keep, a.exclude_rsids, tmp,
                          a.out_tped + ".qc.tsv", a.out_tped + ".vcftools.log")
        os.replace(tmp, a.out_tped)
    finally:
        for f in (tmp, keep):
            if os.path.exists(f):
                os.remove(f)
    print(f"{a.pop} chr{chrom}: {n} samples, {rows:,} sites -> {a.out_tped}")


if __name__ == "__main__":
    main()

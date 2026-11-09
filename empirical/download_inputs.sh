#!/bin/bash
# Download the public empirical inputs into <dir> (default: <repo>/inputs),
# in the layout config/empirical.json expects:
#   <dir>/1kgp/ALL.chr<N>.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz(.tbi)
#   <dir>/1kgp/integrated_call_samples_v3.20130502.ALL.panel
#   <dir>/pyrho_hg19/<POP>/<POP>_recombination_map_hapmap_format_hg19_chr_<N>.txt
#
#   bash empirical/download_inputs.sh [dir] [chromosomes, default "1..22"]
#   bash empirical/download_inputs.sh inputs 22        (chr22 only, for the test config)
#
# Sizes: VCFs ~15 GB for chr1-22 (chr22 ~0.2 GB); pyrho maps 0.9 GB.
# The corpus was built from release v5a, which the NCBI mirror still serves;
# the EBI directory now lists v5b instead. Differences were not assessed.
#
# Two inputs are NOT downloaded here (README section 7):
#   ancestral_dir   chr<N>_ErinG_AA_hg19.tsv.gz  (ARG ancestral alleles, GRCh37)
#   exclude_rsids   exclude_rsid.txt             (dbSNP RsMergeArch merged rsIDs)
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${1:-$REPO/inputs}"; CHR="${2:-$(seq -s ' ' 1 22)}"
KG=https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502
PYRHO=https://zenodo.org/records/11437540/files/hg19_maps.tar.gz
PYRHO_MD5=fc8041c1e427907030359d9a2fd0a7e5
mkdir -p "$DIR/1kgp" "$DIR/pyrho_hg19"

get() {  # url dest
  [ -s "$2" ] && { echo "have $2"; return; }
  echo "get  $1"; wget -q -O "$2.part" "$1" && mv "$2.part" "$2"
}
get "$KG/integrated_call_samples_v3.20130502.ALL.panel" "$DIR/1kgp/integrated_call_samples_v3.20130502.ALL.panel"
for c in $CHR; do
  f=ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  get "$KG/$f" "$DIR/1kgp/$f"; get "$KG/$f.tbi" "$DIR/1kgp/$f.tbi"
done

if [ ! -e "$DIR/pyrho_hg19/.done" ]; then
  get "$PYRHO" "$DIR/hg19_maps.tar.gz"
  got=$(md5sum "$DIR/hg19_maps.tar.gz" | cut -d' ' -f1)
  [ "$got" = "$PYRHO_MD5" ] || { echo "ERROR: hg19_maps.tar.gz md5 $got != $PYRHO_MD5" >&2; exit 1; }
  tar -xzf "$DIR/hg19_maps.tar.gz" -C "$DIR/pyrho_hg19"
  # flatten a leading directory level if the archive has one
  inner=$(find "$DIR/pyrho_hg19" -name 'YRI_recombination_map_hapmap_format_hg19_chr_22.txt' -printf '%h\n' | head -1)
  [ -n "$inner" ] || { echo "ERROR: no YRI chr22 map in the archive" >&2; exit 1; }
  top=$(dirname "$inner")
  [ "$top" = "$DIR/pyrho_hg19" ] || { mv "$top"/* "$DIR/pyrho_hg19/"; rmdir "$top"; }
  rm "$DIR/hg19_maps.tar.gz"; touch "$DIR/pyrho_hg19/.done"
fi
echo "inputs in $DIR; still needed: ancestral_alleles_hg19/ and exclude_rsid.txt (README section 7)"

#!/bin/bash
# Verify the pipeline environment. Exit 0 only if every check passes.
#   bash env/check_env.sh            (activates the env itself)
# Checks: each tool resolves inside the env, links (ldd), and EXITS 0 on a
# real invocation; selscan is v3.0.0 by md5; cosi2 reproduces two production
# replicates byte-for-byte (env/reference_md5.txt).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
source "$HERE/activate.sh"
fail=0
bad() { echo "FAIL: $*"; fail=1; }

for t in coalescent selscan ihsbin xpehhbin python3; do
  p=$(command -v "$t" || true)
  [ -n "$p" ] || { bad "$t not on PATH"; continue; }
  case "$p" in "$DEEPSWEEP_ENV_PREFIX"/*) ;; *) bad "$t resolves outside the env: $p" ;; esac
  miss=$(ldd "$p" 2>&1 | grep -c 'not found' || true)
  if [ "$miss" = 0 ]; then echo "ok  $t -> $p"
  else bad "$t: $miss shared libraries not found ($(ldd "$p" 2>&1 | grep 'not found' | head -1))"; fi
done
for t in selscan ihsbin xpehhbin; do
  "$t" --help >/dev/null 2>&1 || bad "$t --help exits $?"
done
v=$(selscan --help 2>&1 | grep -o 'selscan v[0-9.]*' | head -1 || true)
[ "$v" = "selscan v3.0.0" ] || bad "selscan reports '$v'"
got=$(md5sum "$(command -v selscan)" | cut -d' ' -f1)
[ "$got" = "$SELSCAN_MD5" ] || bad "selscan md5 $got != $SELSCAN_MD5"
python3 -c "import numpy, pandas; assert numpy.__version__ == '2.0.2' and pandas.__version__ == '2.3.1', (numpy.__version__, pandas.__version__)" \
  || bad "python package versions"

# cosi2 determinism against the production binary
W=$(mktemp -d "${TMPDIR:-/tmp}/ds_check.XXXXXX"); trap 'rm -rf "$W"' EXIT
SWEEP='pop_event sweep_mult_standing "sweep" 2 U(0, 2100) E(20) .5 .5-.7 2 U(0, 2100)'
sed "s|^recomb_file .*|recomb_file $REPO/models/test_recom.recom|" "$REPO/models/gravel.par" > "$W/neutral.par"
{ cat "$W/neutral.par"; echo "$SWEEP"; } > "$W/sweep.par"
run_cosi() {  # par seed prefix
  ( cd "$W" && env COSI_NEWSIM=1 COSI_MAXATTEMPTS=5000 coalescent -p "$1" -v -g \
      --genmapRandomRegions --nTriesPerSim 5000 --drop-singletons .25 \
      --tped "$3" -n 1 -M -r "$2" >"$3.log" 2>&1 ) || bad "coalescent exit $? ($(tail -1 "$W/$3.log"))"
}
run_cosi neutral.par 31000000 neut
run_cosi sweep.par 3821012977 sweep
PYTHONPATH="$REPO/pipeline" python3 - "$W" <<'PY'
import os, sys, cms_stats as C
w = sys.argv[1]
for tag, pre in (("neutral_gravel_3mb_seed31000000", "neut"), ("sweep_gravel_3mb_pop2_seed3821012977", "sweep")):
    for p in (1, 2, 3):
        src = os.path.join(w, f"{pre}_0_{p}.tped")
        if os.path.exists(src):
            C.flip_tped(src, os.path.join(w, f"{tag}_pop{p}.tped"))
PY
( cd "$W" && grep -v '^#' "$HERE/reference_md5.txt" | md5sum -c --quiet - ) \
  && echo "ok  cosi2 reproduces the 6 reference TPEDs" || bad "cosi2 output differs from the production binary"

[ "$fail" = 0 ] && echo "ENV OK ($DEEPSWEEP_ENV_PREFIX)" || { echo "ENV CHECK FAILED"; exit 1; }

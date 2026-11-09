#!/bin/bash
# Build the pipeline environment from nothing:
#   1. conda env (env/environment.yml) at $DEEPSWEEP_ENV_PREFIX
#   2. cosi2  (coalescent)        built from source at a pinned commit
#   3. hapbin (ihsbin, xpehhbin)  built from source at a pinned commit
#   4. selscan v3.0.0             the release's prebuilt linux binary, md5-checked
#   5. env/check_env.sh           every tool runs; cosi2 output is byte-checked
#
#   bash env/setup_env.sh                       (env in <repo>/.env)
#   DEEPSWEEP_ENV_PREFIX=/path bash env/setup_env.sh
# Needs: conda, git and internet access. Linux x86-64 (selscan's binary).
# Refuses to touch an existing prefix; remove it yourself to rebuild.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/settings.sh"
P="$DEEPSWEEP_ENV_PREFIX"
JOBS="${JOBS:-$(nproc)}"
[ -e "$P" ] && { echo "ERROR: $P exists; remove it to rebuild" >&2; exit 1; }

command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not on PATH" >&2; exit 1; }
set +u; source "$(conda info --base)/etc/profile.d/conda.sh"; set -u
echo "=== conda env -> $P  ($(date))"
conda env create -p "$P" -f "$HERE/environment.yml"
set +u; conda activate "$P"; set -u
mkdir -p "$P/src"

echo "=== cosi2 @ $COSI2_COMMIT"
git clone -q "$COSI2_REPO" "$P/src/cosi2"
git -C "$P/src/cosi2" checkout -q "$COSI2_COMMIT"
(
  cd "$P/src/cosi2"
  # Flags of the build that made the 2026oct corpus. gcc 14 rejects parts of
  # this 2013-era C++ under -Werror, so every -Werror is removed.
  BASE_CXXFLAGS="-O2 -g -Wno-error=deprecated-declarations -Wno-error=maybe-uninitialized -Wno-error=array-bounds -Wno-error=stringop-overflow -fpermissive -std=c++11"
  chmod +x ./configure
  ./configure CXXFLAGS="$BASE_CXXFLAGS"
  sed -i 's/-Werror[^[:space:]]*//g' Makefile
  for f in -Wno-error=array-bounds -Wno-error=deprecated-declarations -Wno-error=maybe-uninitialized -Wno-error=stringop-overflow -fpermissive -std=c++11; do
    grep -q -- "$f" Makefile || sed -i "/^CXXFLAGS =/s/$/ $f/" Makefile
  done
  make -j"$JOBS"
  install -m 0755 coalescent "$P/bin/coalescent"
)

echo "=== hapbin @ $HAPBIN_COMMIT"
git clone -q "$HAPBIN_REPO" "$P/src/hapbin"
git -C "$P/src/hapbin" checkout -q "$HAPBIN_COMMIT"
(
  mkdir -p "$P/src/hapbin/build" && cd "$P/src/hapbin/build"
  cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_INSTALL_RPATH="$P/lib" \
        -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON ../src
  make -j"$JOBS" ihsbin xpehhbin
  install -m 0755 ihsbin xpehhbin "$P/bin/"
  cp -P libhapbin.so* "$P/lib/"
)

echo "=== selscan $SELSCAN_TAG (prebuilt linux binary)"
git clone -q --depth 1 --branch "$SELSCAN_TAG" "$SELSCAN_REPO" "$P/src/selscan"
got=$(md5sum "$P/src/selscan/bin/linux/selscan" | cut -d' ' -f1)
[ "$got" = "$SELSCAN_MD5" ] || { echo "ERROR: selscan md5 $got != $SELSCAN_MD5" >&2; exit 1; }
install -m 0755 "$P/src/selscan/bin/linux/selscan" "$P/bin/selscan"
rm -rf "$P/src/selscan/.git"

echo "=== check"
bash "$HERE/check_env.sh"
echo "=== done $(date): source env/activate.sh"

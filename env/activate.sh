# Source this (do not execute): activates the pipeline env in a clean shell.
#   source env/activate.sh
# Works under `set -euo pipefail`.
_ds_env_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_ds_env_dir}/settings.sh"
command -v conda >/dev/null 2>&1 || { echo "activate.sh: conda not on PATH" >&2; return 1; }
_ds_had_u=0; case $- in *u*) _ds_had_u=1; set +u ;; esac
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$DEEPSWEEP_ENV_PREFIX" || { echo "activate.sh: cannot activate $DEEPSWEEP_ENV_PREFIX (run env/setup_env.sh)" >&2; return 1; }
[ "$_ds_had_u" = 1 ] && set -u
export HAPBIN_SRC="${DEEPSWEEP_ENV_PREFIX}/src/hapbin"
export DEEPSWEEP_REPO="$(cd "${_ds_env_dir}/.." && pwd)"
unset _ds_env_dir _ds_had_u

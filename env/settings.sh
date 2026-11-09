# Where the conda env is created. Sourced by setup_env.sh and activate.sh.
# Override by exporting DEEPSWEEP_ENV_PREFIX before running either.
_ds_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEEPSWEEP_ENV_PREFIX="${DEEPSWEEP_ENV_PREFIX:-${_ds_repo}/.env}"
unset _ds_repo

# Pinned tool sources (TOOLS.md explains each choice).
COSI2_REPO="https://github.com/broadinstitute/cosi2.git"
COSI2_COMMIT="5da717cb8bf7324fb961e772a0c9042385de381a"   # branch is-251008-drop-imposs-traj
HAPBIN_REPO="https://github.com/evotools/hapbin.git"
HAPBIN_COMMIT="3babf68e0eb4bdfcbda4f56b0821ad4ff141f209"
SELSCAN_REPO="https://github.com/szpiech/selscan.git"
SELSCAN_TAG="v3.0.0"
SELSCAN_MD5="5490bbcd6c7d1c15b3b55aff1f851a26"            # bin/linux/selscan at v3.0.0

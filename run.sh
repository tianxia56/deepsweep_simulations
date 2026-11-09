#!/bin/bash
# One command from a fresh clone to outputs:
#   ./run.sh                              config/config.json (production scale)
#   ./run.sh config/test_100.json         100 sims of each kind, for a laptop
#   ./run.sh config/my_config.json [--workers N]
# Builds the conda env and tools on first use (env/setup_env.sh, ~15-30 min),
# checks them (env/check_env.sh), then runs every stage (run.py).
# Rerunning the same command resumes an interrupted run.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-$REPO/config/config.json}"; shift || true
source "$REPO/env/settings.sh"
if [ ! -x "$DEEPSWEEP_ENV_PREFIX/bin/coalescent" ]; then
  echo "== no environment at $DEEPSWEEP_ENV_PREFIX: building it"
  bash "$REPO/env/setup_env.sh"
fi
bash "$REPO/env/check_env.sh"
source "$REPO/env/activate.sh"
python3 "$REPO/run.py" "$CONFIG" "$@"

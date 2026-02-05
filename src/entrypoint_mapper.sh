#!/usr/bin/env bash
set -euo pipefail

# Always start in repo root
cd /opt/DSE

echo ""
echo "Mapper container ready (conda env: mapper_experiments, Python 3.10)"
echo "Repo: /opt/DSE"
echo ""

# Required defaults for your CLI
CIRCUITS_DIR="${CIRCUITS_DIR:-circuits}"
OUT_FILE="${OUT_FILE:-results/mapper_results.csv}"

# Make sure results dir exists
mkdir -p "$(dirname "$OUT_FILE")"

# Use an array to avoid quoting/spacing bugs
DEFAULT_ARGS=(--circuits-dir "$CIRCUITS_DIR" --out "$OUT_FILE")

# Debug (optional)
conda run -n mapper_experiments python -c "import sys; print('Python:', sys.executable)"
conda run -n mapper_experiments python -c "import pandas as pd; print('pandas OK', pd.__version__)"
echo "Default args: ${DEFAULT_ARGS[*]}"
echo ""

run_default() {
  exec conda run --no-capture-output -n mapper_experiments \
    python -u src/run_mapper_experiments.py "${DEFAULT_ARGS[@]}"
}

# Non-interactive default
if [[ "${RUN_DEFAULT:-}" == "yes" ]]; then
  echo "RUN_DEFAULT=yes -> running default mapper experiment..."
  run_default
fi

read -r -p "Run default mapper experiments now? [y/N]: " ans
if [[ "$ans" =~ ^([yY][eE][sS]|[yY])$ ]]; then
  run_default
else
  echo "Dropping you into a shell with mapper_experiments activated."
  exec bash -lc "source /opt/conda/etc/profile.d/conda.sh && conda activate mapper_experiments && cd /opt/DSE && exec bash -i"
fi

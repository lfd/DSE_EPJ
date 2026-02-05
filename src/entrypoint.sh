#!/usr/bin/env bash
set -e
source /opt/conda/etc/profile.d/conda.sh
conda activate device_experiments

cd /opt/DSE

echo ""
echo "DSE container is ready (conda env: device_experiments)."
echo "Run default experiment now? [y/N]"
read -r ans

if [[ "$ans" =~ ^[Yy]$ ]]; then
  echo "Running default (example): python src/run_device_experiments.py"
  python src/run_device_experiments.py
else
  echo "Dropping you into a shell. You're in /opt/DSE with env activated."
  exec bash --noprofile --norc

fi

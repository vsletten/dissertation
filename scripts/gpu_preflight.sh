#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "$repo_root/qm/quarry/etiquette.py"
nvidia-smi \
  --query-gpu=index,name,memory.total,memory.used \
  --format=csv,noheader

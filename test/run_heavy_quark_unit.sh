#!/bin/bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "${repo_root}/test/bin"

g++ -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I"${repo_root}" \
  "${repo_root}/test/test_heavy_quark_energy_loss.cc" \
  "${repo_root}/HeavyQuarkEnergyLoss.cc" \
  "${repo_root}/Random.cc" \
  -o "${repo_root}/test/bin/test_heavy_quark_energy_loss"

"${repo_root}/test/bin/test_heavy_quark_energy_loss"

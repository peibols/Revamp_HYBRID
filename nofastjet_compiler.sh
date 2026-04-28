#!/bin/bash
set -euo pipefail

# Use the validated PYTHIA 8.315 install used for the main/Moliere closure checks.
PYTHIA_BIN=${PYTHIA_BIN:-/data/yjlee/pythia/pythia8/pythia8315/bin}
PYTHIA_INCLUDE=${PYTHIA_INCLUDE:-/data/yjlee/pythia/pythia8/pythia8315/include}
PYTHIA_LIB=${PYTHIA_LIB:-/data/yjlee/pythia/pythia8/pythia8315/lib}
PYTHIA_SHARE=${PYTHIA_SHARE:-/data/yjlee/pythia/pythia8/pythia8315/share/Pythia8}

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <source-stem>" >&2
  exit 2
fi

g++ -g -std=c++17 -mcmodel=medium -Wall \
    "$1.cc" -o "$1" \
    "${PYTHIA_LIB}/libpythia8.so" -I"${PYTHIA_INCLUDE}" -L"${PYTHIA_LIB}" -Wl,-rpath,"${PYTHIA_LIB}" -lpythia8 -ldl

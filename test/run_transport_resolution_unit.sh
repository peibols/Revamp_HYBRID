#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
binary="${TMPDIR:-/tmp}/mmli_transport_resolution_unit"

g++ -std=c++17 -Wall -Wextra -pedantic -I"${repo_dir}" \
  "${repo_dir}/test/test_transport_resolution.cc" -o "${binary}"
"${binary}"

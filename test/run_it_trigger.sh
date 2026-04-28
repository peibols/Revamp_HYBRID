#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/.." && pwd)
cd "${script_dir}"

njob=1                  # Seed base used across the staged trigger workflow.
Nev=10                  # Number of events.
alpha=0.404             # Strong-coupling energy-loss strength.
kappa=0                 # Gaussian broadening strength.
mode=0                  # 0=strong, 1=radiative, 2=collisional.
tmethod=1               # Tc selector: 145 MeV when set to 1.
cent="0-5"              # Test hydro bundled with this branch.
do_quench="true"        # If false, run the staged trigger path in vacuum mode.

trigger_pt=40.
trigger_eta=1.44
trigger_id=22
sqrts=5020
rpower=2.0

if [[ "${do_quench}" != "true" ]]; then
  alpha=0
  kappa=0
fi

time "${repo_root}/photelsen" "${Nev}" "${njob}" "${sqrts}" "${trigger_pt}" "${trigger_eta}" "${trigger_id}"
time "${repo_root}/plasma_lund" "${njob}" "${alpha}" "${kappa}" "${Nev}" "${mode}" "${tmethod}" "${cent}" "${rpower}" "${njob}" "${sqrts}"
time "${repo_root}/wake" "${njob}" "${njob}" "${Nev}"
time "${repo_root}/hadronizer" "${njob}" "${Nev}" "${njob}"

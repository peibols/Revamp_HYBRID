# Finite-LRES Integration Notes

Branch: `main_moliere_lres_integration`

## Current Milestone

Finite color-resolution energy loss is now available inside the MMI event loop behind:

- `do_lres = true`
- `rpower = <finite value>` or `lres_rpower = <finite value>`

The default path is unchanged:

- `do_lres = false` keeps the existing MMI/Moliere energy-loss path.
- `do_lres = true` with `rpower >= 1e9` falls back to the existing MMI energy-loss path, preserving the already validated `lres -> 0` limit behavior.

The integrated finite-LRES path uses the MMI shower tree, MMI geometry RNG, MMI hydro accessor, MMI in-event wake timing, and MMI Lund hadronization contract. It does not call the old staged executables.

For validation only, the input card also accepts:

- `use_fixed_xy = true`
- `fixed_x = <value>`
- `fixed_y = <value>`

This bypasses Glauber vertex sampling for controlled one-event comparisons. Leave it disabled for physics runs.

## Implementation Summary

- `HYBRID` now reads `do_lres` and `rpower`/`lres_rpower` from the input card.
- `EnergyLoss` dispatches to a finite-LRES resolver when `do_lres = true` and `rpower < 1e9`.
- The resolver builds a vacuum formation timeline, computes daughter resolution times, rewrites effective mother links, then applies the existing MMI energy-loss stepper to resolved objects.
- When `do_elastic = true`, finite LRES owns the propagation timeline and calls Moliere on the currently resolved effective color object. Unresolved daughters do not scatter independently until their resolution time.
- A causal-time guard is applied when a resolved start position would otherwise have `t < |z|`; this prevents invalid proper-time evaluation during the integrated MMI loss-rate step.
- A fixed-vertex validation hook was added in `HYBRID` to isolate vertex RNG differences from finite-LRES physics differences.
- Moliere inverse-CDF kinematics now reject zero/non-finite endpoint draws before evaluating `Cu_`; those endpoints previously produced deterministic `Cu failed 0 0 0 ...` aborts in standalone Moliere as well as in finite-LRES compositions.

## Validation

Pinned PYTHIA setup:

- `/data/yjlee/pythia/pythia8/pythia8315`

Smoke directory:

- `/raid5/data/yjlee/hybrid_dev/test/mmi_lres_integration_smoke_20260430`

Checks:

- Default-path regression, `do_lres = false`: byte-identical to the pre-integration MMI binary for 1 event.
- Default-path regression was repeated after adding the fixed-vertex validation hook: partons and hadrons remain byte-identical to the earlier validated MMI output.
- Current-main core-physics regression, `do_lres = false`, `do_elastic = false`, wake off, 100 events: byte-identical to `origin/main` commit `36151fa` for partons and hadrons.
- Integrated finite-LRES smoke, `do_lres = true`, `rpower = 2.0`, wake off, `kappa = 0.5`: runs successfully for 1 event.
- Integrated finite-LRES reproducibility: two identical 1-event runs are byte-identical.

Current-main 100-event benchmark:

- Directory: `/raid5/data/yjlee/hybrid_dev/test/main_vs_mmli_20260503/run_20260503_124059`
- Reference: clean detached `origin/main` at `36151fa` (`Merge pull request #2 from peibols/main_moliere_integration`).
- Setup: pinned PYTHIA 8.315, averaged hydro, `cent = 0-5`, `kappa = 0.5`, `alpha = 0.404`, `tmethod = 0`, `do_quench = true`, `do_wake = false`, `do_elastic = false`, `do_lres = false`, seed 0, 100 events.
- Result: `pass=True`; return codes OK, hadrons exact, partons exact.
- Output sizes and hashes:
  - Hadrons: 13,688 lines in both branches, SHA256 `4989983e193935a51d3514edc628251a9871f831b33a2dd3cb4c349aaae3073e`.
  - Partons: 1,897 lines in both branches, SHA256 `c28c878dad31bfbbc097d471ef55c76c04f17caa5cb1281ab203ae31ed7d536d`.
- Runtime on the loaded validation node: current main 2.774 s, MMLI 2.719 s.

Finite-LRES + Moliere validation:

- Directory: `/raid5/data/yjlee/hybrid_dev/test/mmli_lres_moliere_integration_20260502`
- Focused endpoint probe after the Moliere kinematic guard: `seed2_probe_20260503_084219`.
- Full matrix after the guard: `moliere_on_run_20260503_084440`.
- Setup: pinned PYTHIA 8.315, averaged hydro, `cent = 0-5`, `do_lres = true`, `lres_rpower = 2.0`, `do_elastic = true`, `compat_moliere_legacy_hydro = true`, 10 seeds, 1 event per seed, two same-seed repeats.
- Matrix cells: wake off/on crossed with broadening off/on, with `kappa = 0` and `kappa = 0.5`.
- Result: all 40 repeat pairs pass; parton and hadron outputs are byte-identical in every cell, with no timeouts.
- Runtime envelope from the passing matrix:
  - wake off, broadening off: average 32.4 s per repeat, maximum 61.6 s
  - wake on, broadening off: average 48.2 s per repeat, maximum 107.6 s
  - wake off, broadening on: average 35.2 s per repeat, maximum 74.3 s
  - wake on, broadening on: average 46.7 s per repeat, maximum 90.9 s

Staged-vs-integrated finite-LRES 1-event matrix:

- Directory: `/raid5/data/yjlee/hybrid_dev/test/mmi_lres_vs_staged_finite_1evt_20260430`
- Cases: wake off/on crossed with broadening off/on, `rpower = 2.0`, `seed_base = 7`.
- Legacy staged finite-LRES versus integrated MMI-LRES does not byte-match in any of the four cells.
- First legacy-mode mismatch is the production vertex: staged legacy samples `X = 5.00696`, `Y = -2.80239`, while MMI samples `X = -3.2331`, `Y = -1.25081`.
- Enabling staged `HYBRID_MMI_COMPAT=1` aligns the vertex, but that staged path calls `runMmiCompatEloss` and bypasses finite-LRES resolution, so it is not a finite-LRES reference.
- Fixed-vertex diagnostic, wake off and `kappa = 0`, confirms the remaining difference is in the finite-LRES loss/timing/output contract rather than the vertex RNG. The same final colored parton count is present, but the rows are not byte-identical.

Historical `untouched_basic` probe:

- Directory: `/raid5/data/yjlee/hybrid_dev/test/main_vs_mmli_20260503/run_untouched_20260503_124402`
- Reference: clean detached `origin/untouched_basic` at `c775edf`.
- Setup: pinned PYTHIA 8.315, legacy `main`, `do_quench = true`, wake compiled on, `kappa = 0.5`, `alpha = 0.404`, `tmethod = 0`, seed 0, 1 event.
- Result: not an exact-reference gate. The legacy branch and MMLI both run, but byte identity fails immediately because the legacy production vertex and medium evolution differ from the MMI-native contract:
  - legacy header: `X = 2.62739`, `Y = 0.759099`
  - MMLI header: `X = -1.69064`, `Y = -2.08526`
- Interpretation: `origin/untouched_basic` is useful historical context, but not the byte-level integration oracle for MMLI. The byte-level production oracle is current `origin/main`/MMI with the MMI-native contract.

Deprecated 100-event wake-on sync attempt:

- Directory: `/raid5/data/yjlee/hybrid_dev/test/mmli_vs_mmi_sync_20260502/run_20260503_093417`
- Setup: MMI versus MMLI, seed 0, 100 events, `do_lres=false` and `lres -> 0`, wake off/on.
- Wake-off branch pairs completed and were byte-aligned.
- The original Python harness exited before writing `summary.tsv`; four wake-on child executables continued as orphaned processes and later exited with logs stopping at event 76/100.
- Result: this run is not used as a validation gate. It was superseded by controlled gates above plus the existing 10-seed wake-on exact sync matrix.

Hashes for the integrated finite-LRES smoke:

- Partons: `f76dac271982f1394729f6a65069dea784156e06109317f5cd76cdade7b91545`
- Hadrons: `5da77894d7381c20a324c8647a233e6a77010c8e44c9acde2ad50322fb96e704`

## Integration Readiness

The integration contract is MMI-native finite LRES:

- LRES controls color-resolution grouping.
- Hydro access, geometry sampling, energy loss, wake, Lund, and output contract remain MMI-native.
- Moliere scattering is applied segment-by-segment to the currently resolved effective color object when `do_elastic = true`.
- Unresolved daughters do not scatter independently until their resolution time; after resolution, daughters propagate and scatter independently.

Closed gates for this contract:

- Default MMI path unchanged when `do_lres = false`.
- `lres -> 0` limit reproduces MMI when `do_lres = true` and `rpower >= 1e9`.
- Current `origin/main` versus MMLI core path is byte-identical for 100 events.
- Finite `rpower = 2.0` with Moliere on repeats byte-identically across the 2x2 wake/broadening matrix at the 10-seed validation scale.
- The Moliere endpoint guard removes deterministic zero-kinematics crashes without changing completed exact-repeat outputs.

Remaining work after merge is physics production, not integration correctness:

- scale selected finite-LRES physics observables to higher statistics;
- profile wake-on runtime tails under production settings;
- compare spectra, energy balance, parton/hadron counts, and wake-sensitive observables for the chosen finite-LRES physics settings.

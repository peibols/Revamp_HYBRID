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

## Status Snapshot: 2026-05-27

This section records the current working understanding of MMLI finite-LRES plus Moliere after the May 27 slide and code-link pass.

### Branch and Documentation State

- Active implementation branch/worktree: `main_moliere_lres_integration` in `/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_integration_clean`.
- The MMLI implementation referenced in the slides is commit `9157d5b2d4b859c6573dd7fd834a54b3134cac45`.
- Overleaf report file: `/raid5/data/yjlee/hybrid_dev/overleaf_69d5314739aed083fc3bbb0a/report/20260527-status.tex`.
- Overleaf has been pushed through commit `4a5ff11 Add LRES tree timeline animation example`.
- A later local edit added code-link buttons to slides 8-12, but that edit has not yet been pushed because the last command-approval/push attempt was blocked by the current usage limit.

### LRES Tree and Timeline Algorithm

- `EnergyLoss::do_lres_eloss_impl` assumes binary shower nodes. Each parton has at most first/second direct daughters for the LRES sibling-pair construction.
- A topology that looks like `1 -> 2 + 3 + 4` must enter the algorithm as sequential binary splittings, for example `1 -> 2 + I`, then `I -> 3 + 4`.
- LRES builds an effective resolution timeline from the shower formation times and the geometric condition for whether the medium can resolve a sibling pair.
- The code enforces ordered effective resolution, so a descendant structure is not allowed to become visible to the medium in a way that is inconsistent with an unresolved ancestor chain.
- The effective object lifetime is then propagated segment-by-segment. Before resolution, the medium sees one effective color object. After resolution, daughters propagate independently.

### Moliere with LRES: Implemented Modes

Mode precedence in the current implementation is:

1. Mode D: `do_Moliere_dynamic_daughter_unresolved_resolution = true`
2. Mode C: `do_Moliere_dynamic_unresolved_resolution = true`
3. Mode B: `do_Moliere_on_unresolved_partons = true`
4. Mode A: default coherent unresolved parent

Mode A:

- During an unresolved finite-LRES interval, Moliere scattering is applied to the coherent effective parent object.
- One accepted coherent scattering creates one recoil/hole source.
- This is the backward-compatible default when all newer unresolved-Moliere options are false.

Mode B:

- During the same unresolved interval, the current parent energy is projected onto the two vacuum daughter directions.
- Each unresolved daughter is propagated separately through Moliere.
- At the end of the unresolved interval, the daughters are recombined into the effective parent four-vector expected by the surrounding LRES algorithm.

Mode C:

- Candidate scatterings are generated from the coherent parent.
- Each parent-level candidate is tested with `q_perp * d_perp > c_res`.
- If the test fails, the kick remains coherent on the parent.
- If the test passes, the pair is marked elastically decohered and the parent-level kick is assigned to one daughter probabilistically, weighted by daughter energy.
- This is implemented, but the struck daughter assignment is approximate because the candidate was generated from the parent, not from a particular daughter.

Mode D:

- Candidate scatterings are probed separately from each daughter using copied `numrand` states and `propagate_segment_with_scattering_callback` with `StopBeforeApply`.
- The earliest daughter candidate is selected.
- Before testing that candidate, the coherent parent is propagated from the current time to the candidate time with the normal HYBRID `loss_rate`, preserving continuous unresolved energy loss.
- The candidate is tested with `q_perp * d_perp > c_res`.
- If the test fails, the daughter-probe sampled kick is applied coherently to the parent and exactly one recoil/hole source is created. This is an implemented approximation: the candidate came from a daughter probe, but its wavelength is treated as too long to resolve the dipole.
- If the test passes, the struck daughter receives the kick, the pair is marked elastically decohered, both daughters are propagated independently for the rest of the original unresolved interval, and the daughters are recombined at the segment boundary.
- The daughter probes are side-effect-free for momenta/recoilers/holes, but they intentionally consume/advance RNG through copied `numrand` objects and then update `nr_`.

### Dynamic Resolution Helpers

- The dipole size helper currently uses an approximate transverse separation:
  `d_perp ~= |(v_perp,1 - v_perp,2) * (t_candidate - t_split)|`.
- The resolution condition is `q_perp * d_perp > c_res`.
- The runtime parameter `moliere_unresolved_resolution_c` defaults to `1.0`.
- The dynamic-resolution test affects only the unresolved elastic/Moliere handling. It does not change hydro, the LRES timeline construction, the standard energy-loss formula, or downstream analysis selection.

### Event Display and Tree/Timeline Dumps

- Existing HYBRID diagnostic options:
  - `dump_hybrid_evolution_history = true` writes time-ordered LRES/Moliere/response records to a TSV file.
  - `doEventDisplay = true` writes ROOT TTrees for propagated parton time segments and detailed records.
- A standalone plotting utility was added locally at:
  `/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_integration_clean/test/plot_lres_tree_timeline.py`.
- That utility converts a HYBRID history TSV into an independent JSON tree/timeline dump plus PNG/GIF animation.
- Example generated files are stored at:
  `/raid5/data/yjlee/hybrid_dev/test/lres_tree_timeline_example/`.
- The corresponding Overleaf assets were pushed under:
  `eventDisplay/lres_tree_timeline_example/`.

### Slides

- `20260527-status.tex` now documents the LRES tree/timeline construction, effective objects, propagation handoff, energy-loss stepper, and Moliere modes A-D.
- Slides 2-6 already had working GitHub code buttons after switching to `\beamergotobutton{code}` links.
- Slides 8-12 have local code-link buttons added for the mode A-D algorithms and implementation contract. These were compiled successfully locally, but remain unpushed as of this snapshot.
- The local compile produced `20260527-status.pdf` successfully with 14 pages and only a small overfull warning on the Mode-D outcome slide.

### Open Caveats

- Mode D is closer to a daughter-candidate interpretation than Mode C, but it is still not a full analytic merged Poisson process for unresolved dipoles.
- Failed Mode-D tests coherently apply a kick that was sampled from a daughter probe. This is deliberate in the current implementation and should be described as an approximation.
- Mode E follow-up implemented locally: when a coherent object opens, daughter momenta are now materialized by rotating the vacuum daughter directions into the live parent axis and conserving the live parent energy. If no elastic candidate decoheres the parent before the normal finite-LRES boundary, the boundary split now seeds live daughter states so coherent parent deflections are inherited by later descendants. The residual spatial-momentum mismatch from opening an on-shell coherent parent is an explicit approximation to validate.
- Mode E follow-up implemented locally: recursive `q_perp d_perp` tests now use live projected daughter positions when available (`qperp_dperp_live` in the history output) and fall back to the old vacuum estimate only if projection fails (`qperp_dperp_vac_fallback`). Diagnostics count live tests and fallbacks.
- Remaining TODO for Mode E: the sampled momentum-transfer delta is still taken from the frontier probe and then mapped upward if the dipole test fails. This keeps one kick and one recoil/hole source, but it is not yet a fully analytic merged Poisson process for the coherent object.
- Remaining TODO for Mode E: audit whether deterministic frontier order (`4, 5, 3`, for example) biases candidate selection or RNG consumption. The robust target is an order-independent candidate scheduler, not a cosmetic shuffle.
- TODO validation study: for the nested example `1 -> 2 + 3`, `2 -> 4 + 5` (final frontier `4,5,3`), compare how the angular distribution relative to the original parent-1 direction changes with and without color-coherence treatment. Track angles such as `DeltaR(4,1)`, `DeltaR(5,1)`, `DeltaR(3,1)`, and the effective-subtree axes before/after coherent parent kicks, then compare coherent propagation, independent daughter propagation, and recursive Mode-E-style coherence.
- The event-display tree/timeline animation is currently a diagnostic visualization, not a physics validation observable.
- The slide source has local changes not yet pushed to Overleaf after the slide 8-12 code-link update.

## Status Snapshot: 2026-06-06

This section records the latest slide/documentation state after preparing the June status deck.

### Overleaf Slide Deck

- The June status deck has been renamed from `report/20260606-status.tex` to `report/20260609-status.tex`.
- The displayed slide date is now `June 9, 2026`.
- The rename/date update was pushed to Overleaf in commit `600c574 Rename status deck to 20260609`.
- The deck compiles locally as `report/20260609-status.pdf` with 19 pages. The generated PDF and aux/log files are local build products and were not committed.

### Current Slide Content

- The deck retains the finite-LRES implementation slides: tree/timeline construction, effective objects, propagation handoff, energy-loss stepper, and current Moliere modes A-D.
- The mode A-D slides include clickable GitHub code buttons pointing to the MMLI implementation commit used for the documentation.
- Two Mode-E planning slides were added and pushed in Overleaf commit `322d432 Add Mode E implementation path slides`.

### Mode E Direction

Mode E is intended to be the fully Korinna-like dynamic colour-coherence implementation, not just another per-segment unresolved-Moliere option.

The required changes are:

1. Promote colour coherence to explicit event state: active dipoles/coherent systems, members, formation times, positions, and decoherence flags.
2. Replace the current per-segment A-D branch with an interaction-level candidate scheduler.
3. Test each actual sampled elastic/Moliere momentum transfer against the relevant dipole size, using `q_perp * d_perp > c_res`.
4. Failed test: apply one coherent kick to the coherent system and create one recoil/hole source.
5. Passed test: apply the kick to the struck parton, create one recoil/hole source, and update the coherence graph.
6. After decoherence, subsequent elastic candidates in the same original geometric LRES interval must use the updated active resolved objects.
7. A minimal MMLI Mode E can record decoherence while keeping the PYTHIA vacuum shower fixed, but a fully Korinna-like implementation also needs the shower-evolution consequence of decoherence: an in-medium shower hook, emission veto/reweighting, or regenerated constrained shower.

### Important Caveat

Mode D is the closest implemented approximation today. It uses daughter-level probes and selects the earliest candidate, but a failed dynamic test still applies a daughter-probe sampled kick coherently to the parent. This is a deliberate approximation and should not be described as the full Korinna method.

## Status Snapshot: 2026-06-29

This section records the Mode-E follow-up after the June 16 Dani/Krishna discussion and the June 29 local implementation pass.

### Mode E Changes Implemented Locally

Mode E is the recursive unresolved-Moliere mode enabled by:

- `do_Moliere_recursive_unresolved_resolution = true`

The current local implementation now includes the following follow-up changes beyond the first Mode-E commit:

1. Live daughter materialization from the coherent parent:
   - When a coherent object opens, daughter momenta are constructed by rotating the vacuum daughter directions into the current live parent axis.
   - The live parent energy is partitioned by the vacuum daughter energy fractions.
   - The normal finite-LRES boundary also seeds live daughter states if no elastic candidate decoheres the parent before the boundary, so coherent parent deflections are inherited by later descendants.

2. Live projected `d_perp` in recursive tests:
   - Recursive `q_perp d_perp` tests first project both siblings from the current active tree state to the candidate scattering time.
   - History records use `qperp_dperp_live` when this live projection succeeds.
   - The older vacuum estimate is retained only as a fallback and is labeled `qperp_dperp_vac_fallback`.
   - Diagnostics count `n_recursive_live_dperp_tests` and `n_recursive_vacuum_dperp_fallbacks`.

3. Branch-local frontier probe RNGs:
   - The previous Mode-E candidate scan used one copied `numrand` stream in deterministic tree order. For a frontier such as `4, 5, 3`, the first branch got the first random draws, the second branch got the next random draws, etc.; this made the sampled candidates depend on an arbitrary traversal order.
   - The current local code now assigns each frontier probe a deterministic branch-local seed keyed by the base HYBRID seed, event id, LRES segment id, Mode-E iteration id, active ancestor, and probe id.
   - The Moliere elastic sampler also had to be scoped explicitly, because `gen_particles` uses the global `std::default_random_engine generator` from `Distributions.hpp`, not `numrand`. Mode E probes now snapshot, seed, and restore that elastic generator around each branch-local probe.
   - Frontier ids are sorted before probing, and equal-time candidate ties are broken by probe id.
   - The history dump records each selected probe candidate as `recursive_probe_candidate` with a `modeE_branch_local_seed=...` note.
   - New diagnostics count `n_recursive_frontier_probe_batches`, `n_recursive_frontier_probe_objects`, `n_recursive_frontier_permutation_checks`, and `n_recursive_frontier_permutation_mismatches`.

4. Coherent upward application location:
   - If the recursive tests map a frontier-probe candidate upward to an unresolved coherent object, the code applies one kick and one recoil/hole source at the live coherent object/subtree state.
   - This preserves the intended one-source behavior for unresolved coherent scatterings.

### Important Approximations That Remain

- The momentum-transfer delta is still sampled from a frontier probe. If the `q_perp d_perp` tests fail, the same sampled kick is mapped upward to the coherent parent/object. This is a controlled implementation approximation, not yet an analytic merged coherent-object Poisson process.
- The branch-local RNG scheduler now removes dependence on arbitrary frontier traversal order for Mode-E probes, including the Moliere elastic `Distributions.hpp` generator. It is still a stochastic convention rather than a derivation of the exact coherent-system elastic rate.
- A coherent massless/on-shell parent cannot generally be opened into two separated massless/on-shell daughters while preserving the exact full four-vector, live opening angle, and causal daughter kinematics simultaneously. The current implementation conserves the live parent energy and direction/deflection, then treats the residual spatial-momentum mismatch as a validation item.
- The PYTHIA vacuum shower is still fixed. A fully Korinna/JEWEL-like treatment would also need an in-medium shower hook, emission veto/reweighting, or constrained regeneration after elastic decoherence.

### Validation Run After Branch-Local Probe Update

Build:

- Compiled `main` successfully with pinned PYTHIA 8.315 from `/raid5/data/yjlee/hybrid_dev/test/tmp_pythia8315_validate_mmi_revert/pythia8315`.
- Existing compiler warnings only; no new compile errors.

Smoke checks under `/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_integration_clean/test/tmp_modeE_recursive_smoke/`:

- `modeE_smoke20_branch_rng.input`, `moliere_unresolved_resolution_c = 1.0`, 20 events: completed successfully.
  - `n_unresolved_segments_dynamic = 130`
  - `n_unresolved_candidate_scatters = 8`
  - `n_unresolved_coherent_scatters = 8`
  - `n_unresolved_resolving_scatters = 0`
  - `n_unresolved_pairs_elastically_decohered = 0`
  - `n_recursive_frontier_probe_batches = 137`
  - `n_recursive_frontier_probe_objects = 183`
  - `n_recursive_coherent_applications = 8`
  - `n_recursive_tree_updates = 130`
  - `n_recursive_live_dperp_tests = 0`
  - `n_recursive_vacuum_dperp_fallbacks = 0`
- `modeE_smoke20_c0_branch_rng.input`, `moliere_unresolved_resolution_c = 0.0`, 20 events: completed successfully with the same aggregate scheduler diagnostics in this sample.

Interpretation:

- These two smoke samples validate the branch-local frontier candidate scheduler, history output, and coherent upward-application path.
- They did not find a candidate requiring a nested sibling `q_perp d_perp` test, so they do not yet validate a resolving Mode-E nested-dipole event.
- A targeted event search is still needed for the Dani/Krishna topology where a candidate appears on `4`, `5`, or `3` while the ancestor object is still unresolved.

Additional permutation diagnostic after finding the hidden global-generator issue:

- A first permutation check deliberately reran each already-collected frontier in reversed order. It found `3` mismatches in `38` checks, showing that branch-local `numrand` alone did not control Moliere elastic candidate generation.
- Root cause: `gen_particles` samples from the global `std::default_random_engine generator` in `Distributions.hpp`.
- Fix: Mode E probes now snapshot, seed, and restore that elastic generator for each branch-local probe.
- `modeE_smoke20_permutation_fix.input`, `moliere_unresolved_resolution_c = 1.0`, 20 events: completed successfully.
  - `n_unresolved_segments_dynamic = 129`
  - `n_unresolved_candidate_scatters = 5`
  - `n_unresolved_coherent_scatters = 5`
  - `n_unresolved_resolving_scatters = 0`
  - `n_recursive_frontier_probe_batches = 133`
  - `n_recursive_frontier_probe_objects = 178`
  - `n_recursive_frontier_permutation_checks = 32`
  - `n_recursive_frontier_permutation_mismatches = 0`
  - `n_recursive_tree_updates = 129`

### Remaining TODOs

1. Find or construct a Mode-E event with a true nested-dipole candidate:
   - Example topology: `1 -> 2 + 3`, followed by `2 -> 4 + 5`, with a candidate on `4`, `5`, or `3` while the ancestor object is still unresolved.
   - Run A/B/C/D/E side-by-side and record which object receives the kick, which recoil/hole source is created, and how the active coherent groups change.

2. Validate the angular effect of coherence:
   - Track `DeltaR(4,1)`, `DeltaR(5,1)`, `DeltaR(3,1)`, and effective-subtree axes before/after coherent parent kicks.
   - Compare coherent propagation, independent daughter propagation, and recursive Mode-E propagation.

3. Quantify branch-local scheduler stability:
   - Done for the current implementation: the diagnostic reversed-frontier permutation check is implemented and gives `32` checks with `0` mismatches in the latest 20-event smoke.
   - Still future work: compare branch-local Mode E statistically against any future analytic merged coherent-object scheduler.

4. Quantify four-momentum closure at LRES openings:
   - Record the residual between the live coherent parent four-vector and the sum of materialized daughter four-vectors.
   - Decide whether the current energy/direction-preserving split is sufficient for diagnostics, or whether a different off-shell bookkeeping object is needed.

5. Define the full Korinna-like path beyond current MMLI:
   - Current Mode E records elastic decoherence while keeping the PYTHIA vacuum shower fixed.
   - A complete treatment needs an explicit interface for shower evolution after decoherence.

# OO Plan B: arXiv:2509.19430v2

This directory contains the reproducible pre-equilibrium-medium construction
used for the paired O16+O16 HYBRID campaign.

The paired local-sample ROOT schema, wake-label mapping, anti-kT/Soft-Drop
definitions, and conversion command are documented in
[`ROOT_TREE_SCHEMA.md`](ROOT_TREE_SCHEMA.md).

## Published Prescription

The implementation follows arXiv:2509.19430v2, Eqs. (2)-(4):

- solve the implicit effective-temperature equation using the QCD kinetic
  energy attractor;
- anchor to the first hydro slice with the published viscous correction and
  `eta/s=0.12`;
- use the three-flavor conformal equation of state;
- scale transverse flow linearly with `tau/tau_hyd` and use Bjorken
  longitudinal flow;
- apply no quenching before `tau_min`.

The required attractor table is
`reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv`. Its SHA256 is:

```text
1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0
```

It is reconstructed from the public QCD EKT data associated with
arXiv:1908.02866, DOI `10.4119/unibi/2939684`. The production runner rejects a
different checksum, a missing table, `eta/s != 0.12`, or a disabled viscous
anchor.

## Validation

From this directory, with a staged hydro tree available at `staged_hydro/`:

```bash
python3 analysis/test_2509_prehydro_public.py -v
python3 analysis/test_analyze_oo_prehydro_pair.py -v
python3 analysis/test_convert_oo_paired_to_root.py -v
python3 analysis/test_plot_oo_jet_variables.py -v
python3 analysis/test_plot_oo_jet_charge_response.py -v
python3 analysis/test_plot_oo_jet_paired_substructure.py -v
python3 analysis/test_plot_oo_jet_formation_time.py -v
python3 analysis/test_summarize_oo_jet_pt_slices.py -v
python3 cern_support/test_monitor_oo_1m_prehydro_raa.py -v
python3 cern_support/test_run_chunk_job.py -v
python3 cern_support/test_supervise_oo_prehydro_only.py -v
python3 cern_support/test_supervise_oo_50k.py -v
```

The test covers the attractor checksum and known values, implicit-equation
residual, interpolation accuracy, onset, flow, conformal EOS, and pinned
production options.

`cern_support/monitor_oo_1m_prehydro_raa.py` supports an AA campaign with a
separate retained pp snapshot through `--pp-local-eos`. Milestone analysis
always enables strict one-event paired-AA completeness. If more than one
milestone has been crossed since the preceding poll, it analyzes only the
highest one rather than rescanning the same snapshot.

The RAA analyzer reconstructs each PYTHIA run's final `sigmaGen` and
`weightSum` from `HYBRID_Hadrons.out`. Independent one-event AA runs are
combined with the PYTHIA 8.315 `PythiaParallel` convention,
`sigmaGen = sum(weightSum_r * sigmaGen_r) / sum(weightSum_r)`, before the
merged weighted histogram is normalized once. Applying `weight * sigmaGen`
event by event is invalid for this job layout and produces artificial
high-pT suppression.

The jet-variable analyzer accepts an optional inclusive upper boundary through
`--pt-max`. The three campaign slices use the disjoint convention
`30 < pT <= 50`, `50 < pT <= 80`, and `pT > 80 GeV`. Run the analyzer once per
slice, then pass the three metadata files and the inclusive `pT > 30` metadata
to `analysis/summarize_oo_jet_pt_slices.py`. The summary refuses gaps,
overlaps, inconsistent normalization metadata, or sliced cross sections that
do not close to the inclusive result.

The `Zg` and `Rg` panels reserve their first, shaded bin for jets with
`SoftDropValid=0`. The sentinel bin has the same width as one physical bin, so
its plotted integral is the failed-Soft-Drop jet cross section. The TSV labels
it `softdrop_failed`; metadata and summary tables store weighted pass and fail
fractions. Physical `Zg`/`Rg` moments continue to use successful Soft Drop jets
only.

All jet-variable shape panels use the per-variant normalization
`(1/sigmaJetSelected) dSigma/dx`. Their lower panels are ratios of those
normalized shapes, not absolute-yield ratios; the selected-jet denominator is
recomputed inside each paired delete-one-event jackknife replica. Absolute
cross sections and integrated jet yields remain available in the TSV/JSON
audit and are used for jet RAA.

`analysis/plot_oo_jet_charge_response.py` tests whether the incremental
prehydro shift depends on a single-core versus many-core fragmentation proxy.
Its primary proxy is the wake-excluded effective multiplicity
`N_eff=(sum pT)^2/sum(pT^2)` of normal hadrons. It selects and bins using the
no-prehydro jet, then measures `1-pT(Plan B)/pT(no prehydro)` for the matched
jet. Quark/gluon subsets require the same outgoing hard-marker PDG ID in both
variants. This is a final-state correlation study; a causal count of charges
active at hydro start requires retaining the shower formation timeline.

`analysis/plot_oo_jet_paired_substructure.py` performs the migration-safe
follow-up for Soft Drop and momentum dispersion. It selects only the
no-prehydro jet in pT and eta, follows its one-to-one pair-axis match without
placing a second pT cut on the Plan-B jet, and reports the full fail/pass
transition matrix. It also measures paired changes in `PtD`, normal-only
`PtD`, positive and signed multiplicity, both-pass `Zg` and `Rg`, girth,
maximum `kT`, and the R=0.4/0.8 all-split and hardest-split formation-time
estimators. Their ratio-of-means shifts use the common matched set. It also
reports
`1-pT(Plan B)/pT(no prehydro)` for all, quark-tagged, gluon-tagged,
Soft-Drop, and normal-effective-multiplicity subsets. The biased-PYTHIA event
weight is applied once. Delete-one-hard-event jackknife errors are nominal;
delete-one-`hydroIndex` block errors and event-weight concentration diagnostics
are retained in TSV and JSON outputs.

The single-like and many-like classes are defined from no-prehydro
`NormalPtD` through `NormalEffectiveMultiplicity=1/NormalPtD^2`. Their paired
`PtD` shifts are therefore migration/closure diagnostics with an intrinsic
category-boundary correlation, not independent evidence for a charge-count
effect. The incremental matched-pT loss versus that fixed classifier is the
physics-facing comparison.

RAA plots use a logarithmic pT axis and the fixed edges
`4,5,7,10,14,24,36,50,80,150 GeV`. These edges were selected once from the
8426-pair pilot by requiring the worse relative statistical uncertainty of the
two variants to stay near 25%; the observed range was 18.7--24.9%. The edges
remain fixed for later milestones rather than being retuned on each snapshot.

A V2 continuation uses globally unique task IDs and hard seeds. Generate its
manifest with `cern_support/extend_oo_v2_task_manifest.py`; for the second 50k,
the ranges are task IDs `50000--99999`, hard seeds `950000--999999`, and
milestone blocks `11--20`. The hydro assignment is repeated exactly, preserving
the first campaign's Ncoll-weighted exposure to all 500 hydro events. The
analyzer and supervisor accept shifted half-open task ranges through
`--aa-task-start` and `--aa-task-limit`.

After both 50k supervisors write their strict-completion markers, run
`analysis/run_oo_v2_100k_final_analysis.sh`. It combines the two local EOS
snapshots against `aa_task_manifest_combined_100k.tsv`, rejects duplicate or
missing task IDs, requires exact 100,000-row manifest closure, and produces one
100k ROOT file. Hadron RAA, jet RAA for R=0.1/0.2/0.4/0.8, inclusive and
pT-sliced jet spectra/substructure, charge-response, and migration-safe paired
substructure outputs are generated only after those gates pass.
`analysis/watch_oo_v2_100k_final_analysis.sh` is the restart-safe unattended
entry point: it requires both markers to report `accepted_pairs=50000`, exits
immediately if a final `status=PASS` marker already exists, and otherwise runs
the same combined analysis when both halves are ready.

## Prehydro-Only Alpha Retuning Leg

`run_oo_validation_chunk.py --run-prehydro-only` runs exactly one Plan-B AA
variant. It is mutually exclusive with `--run-prehydro-pair`; no no-prehydro
baseline is generated. The output keeps the conventional `with_prehydro`
variant label and records the requested alpha in both `summary.tsv` and
`task_NNNNN_prehydro_only_summary.tsv`.

The 100k alpha=0.335 campaign is prepared and submitted with:

```bash
production/oo_planb_2509/cern_support/submit_oo_v2_prehydro_only_100k.sh
```

The campaign builder requires the combined task manifest SHA256
`2317bd840ee8fffcc71f4b216a34ec4c22d5d7fa91eddb78dc81c1cd7f0a9a5e`.
This fixes task IDs `0--99999`, hard seeds `900000--999999`, and the same
hydro assignment used by the existing two-leg V2 sample. It repackages the
existing V2 executable without rebuilding it and requires executable SHA256
`0d4c68b6ded87379e598e41670dd0bf8f8a39c234f0a9176b1b73ce6e2d88649`.
Only the Python runner and embedded 100k manifest change.

Ten submit files with 10,000 jobs each respect the CERN per-submission limit.
Condor output, error, and event logs are disabled because each wrapper uploads
its status and physics archive directly to EOS. The jobs and their strict
retries use priority `100`, above the V2 reference retries at priority `0`.
This changes ordering for newly available slots but does not preempt reference
jobs that are already running. Run
`cern_support/supervise_oo_prehydro_only.py` against the campaign work area to
audit every archive and resubmit only missing, failed, or malformed task IDs.
The strict audit rejects a wrong alpha, wrong seed or hydro provenance, a
nonzero return code, or any archive containing an accidental baseline leg.

For RAA, jet spectra, and substructure comparisons, join the alpha=0.335 leg
to the existing no-prehydro and alpha=0.37 archives by `task_id`. Before using
a joined task, require identical hard seed, hydro slot/event/Ncoll, hydro
payload checksum, PYTHIA event weight, and hard marker. The existing one-million
event pp denominator is reused; this campaign submits no pp jobs.

`analysis/freeze_oo_v3_snapshot.py` freezes the strict intersection of one or
more disjoint paired-reference mirrors and the alpha=0.335 mirror. It records
the exact task list and hashes, then hard-links only the common status/archive
pairs into separate reference and alpha trees. Feed that frozen task list to
`analysis/convert_oo_v3_to_roots.py`; it builds two pair-aligned ROOT files,
`no/alpha=0.37` and `no/alpha=0.335`, and aborts on any disagreement in task,
hydro provenance, event metadata, hard vertex, or outgoing hard-parton
markers. After running the standard jet-RAA, jet-variable, and paired-
substructure analyzers on both ROOT files, `analysis/make_oo_v3_comparison.py`
produces the three-way hadron/jet overlays and matched substructure summary.
Completion-order snapshots must retain the provisional diagnostic label until
the final common 100k task set is available.

Run the complete V2-style three-way jet-spectrum suite with:

```bash
analysis/run_oo_v3_jet_spectra.sh \
  ROOT_NO_ALPHA037.root ROOT_NO_ALPHA0335.root OUTPUT_DIR EXPECTED_EVENTS
```

The wrapper analyzes the inclusive selections `pT > 20` and `pT > 30 GeV`
and the disjoint selections `20--30`, `30--50`, `50--80`, and `>80 GeV` for
anti-kT radii 0.1, 0.2, 0.4, and 0.8. It produces three-way kinematics and
substructure overlays for every selection and radius (48 detailed plots) plus
an integrated momentum-slice summary. The merge requires the no-prehydro
histograms from the two aligned ROOT pairs to agree exactly. The summary also
requires the four disjoint intervals to close to the inclusive `pT > 20 GeV`
cross section for every radius and all three variants.

The schema-v6 ROOT files retain the schema-v5 formation-time branches and add
the signed `TotalMult = NNormal + NPositiveWake - NNegativeWake` jet branch.
Run the formation-time analysis separately on each aligned ROOT
pair:

```bash
python3 analysis/plot_oo_jet_formation_time.py \
  --input-root ROOT_NO_ALPHA037.root \
  --out-dir OUTPUT_DIR/alpha037 \
  --prefix oo5360_alpha037_formation_time

python3 analysis/plot_oo_jet_formation_time.py \
  --input-root ROOT_NO_ALPHA0335.root \
  --out-dir OUTPUT_DIR/alpha0335 \
  --prefix oo5360_alpha0335_formation_time
```

Pass both spectra, summary, correlation, and metadata outputs to
`analysis/plot_oo_v3_jet_formation_time.py` for the three-way merge. The merge
requires exact agreement of the duplicated no-prehydro leg. Results cover
R=0.4 and R=0.8 in `(30,50]`, `(50,80]`, and `(80,infinity)` GeV. The C/A
tree contains normal plus positive-wake hadrons; negative wake and hadronized
holes affect the corrected jet selection through jet-level 4MomSub but do not
define a negative-subtracted nonlinear tree. Both all-tree and global
hardest-`kT` estimators use paired delete-one-run jackknife errors. This is a
final-constituent formation-time estimator, not generator shower history.

For a clearly labeled completion-order provisional snapshot before exact
closure, pass each disjoint local EOS mirror as a repeated `--source` to
`analysis/freeze_oo_v2_snapshot.py`, together with the combined manifest. The
freeze rejects overlapping task IDs, records per-source file counts and hashes,
and hard-links only strict paired archives into one immutable `local_eos` tree.
Run `analysis/run_oo_v2_provisional_snapshot.sh` on that frozen directory; its
marker remains `PROVISIONAL_DIAGNOSTIC_NOT_UNBIASED` by construction.

The wrapper separates physics output storage from immutable input storage.
`EOS_BASE` receives the new status and output archives, `PAYLOAD_EOS_BASE`
provides the pinned PYTHIA 8.315 and hydro archives, and
`RUNTIME_PAYLOAD_EOS_BASE` can provide a continuation-specific runtime with its
embedded manifest. This avoids copying the roughly 1 GB hydro payload set for
each continuation campaign.

For monitoring more than one live AA root, repeat
`--sync-additional-aa EOS_BASE LOCAL_EOS` together with the matching
`--additional-aa-local-eos LOCAL_EOS`. The monitor synchronizes status first
and downloads outputs only when a milestone is eligible. Synchronization uses
incremental rsync, so unchanged archives are retained locally while repaired
status and output files are refreshed. Static retained roots remain ordinary
`--additional-aa-local-eos` inputs and are never overwritten.

Large Condor arrays default to `TOLERATE_CHUNK_FAILURES=true`. A failed process
still uploads `status=failed` and its original exit code, but the wrapper exits
successfully to Condor so DAGMan does not remove unrelated processes. Milestone
analysis accepts only paired `status=success` archives; failed or missing task
IDs are resubmitted explicitly after the array drains. The campaign-specific
`cern_support/supervise_oo_50k.py` waits for both guarded roots to drain,
incrementally downloads them, audits their expected task IDs, runs the strict
content analyzer, and resubmits missing, failed, or malformed live-root IDs
until 50,000 parsed pairs are local.

The V2 supervisor removes terminal held records after auditing EOS, including
records whose physics archive is incomplete and must be retried. Retry submit
files set `transfer_output_files = ""` because the wrapper uploads the guarded
status and physics archive directly to EOS; this avoids a redundant Condor
transfer into AFS after successful jobs and prevents an AFS quota failure from
stalling final retries.

## Reproducibility Boundary

This reproduces every pre-equilibrium-medium operation stated publicly in the
paper. It does not reproduce the authors' complete semi-analytic energy-loss
calculation, and the paper does not publish its exact runtime attractor object.
An author reference output is required before claiming author-code or bitwise
identity.

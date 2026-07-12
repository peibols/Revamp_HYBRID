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
python3 analysis/test_summarize_oo_jet_pt_slices.py -v
python3 cern_support/test_monitor_oo_1m_prehydro_raa.py -v
python3 cern_support/test_run_chunk_job.py -v
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
`PtD`, their ratio-of-means shifts on the common matched set, and
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

# OO Paired ROOT Trees

This document defines the ROOT conversion for the paired no-prehydro and
Plan-B-prehydro O16+O16 HYBRID samples. The converter reads the local CERN
archives directly; it never extracts a second copy of the event sample.

The implementation is:

```text
analysis/convert_oo_paired_to_root.py
analysis/oo_root_tree_writer.cc
```

The Python front end validates archive completeness and streams accepted pairs
to the C++ writer. The C++ writer uses ROOT and FastJet. Only archives with a
local `status=success`, exactly one event in both variants, successful pair
summary rows, and matching event number, PYTHIA weight, `sigmaGen`, hard
vertex, seed, and hydro index enter the ROOT file. Every rejection is recorded
in the audit TSV. For v2, pass `--aa-task-manifest`; the converter then also
requires the seed, hydro slot, hydro event ID, Ncoll, and hydro-payload SHA256
in both variant summaries to match the immutable task manifest.

## File layout

One ROOT file contains the complete paired sample:

| Object | Entry definition |
| --- | --- |
| `noPrehydro/Hadrons` | one entry per accepted no-prehydro event |
| `withPrehydro/Hadrons` | one entry per accepted prehydro event |
| `noPrehydro/Jets` | one entry per accepted no-prehydro event |
| `withPrehydro/Jets` | one entry per accepted prehydro event |
| `Pairs` | one entry per strict no/with-prehydro pair |
| `metadata/Sources` | mapping from `sourceIndex` to the local input root |
| `metadata/Conversion` | total pair, hadron, and jet counts |

All four event trees have identical entry ordering. `pairId` is therefore the
preferred join key, not the TTree entry number. It is encoded as
`(sourceIndex << 48) | chunkId`, which keeps overlapping chunk IDs from
different CERN campaigns distinct.

Common event branches include `pairId`, `sourceIndex`, `chunkId`,
`hydroIndex`, `seed`, `eventNumber`, `prehydroEnabled`, `eventWeight`,
`sigmaGen`, `hardX`, and `hardY`. `eventWeight` is the biased-PYTHIA event
weight. It must be combined with the run-level `sigmaGen/weightSum`
normalization used by `analyze_oo_prehydro_pair.py`; it must not be interpreted
as a complete per-event cross-section weight.

## Hadron trees

The requested branches are jagged vectors:

```text
hadronPt  hadronEta  hadronPhi  hadronStatus  eventWeight  hadronID
```

The conversion also stores `hadronRawLabel`, `hadronMass`, `hadronE`,
`hadronPx`, `hadronPy`, and `hadronPz`, so the original four-vector and the two
different kinds of negative contribution are recoverable.

The status convention is:

| `hadronStatus` | HYBRID raw label | Meaning |
| ---: | ---: | --- |
| `0` | `0` | normal quenched hadron |
| `+1` | `1` | positive wake hadron |
| `-1` | `2` | negative wake hadron |
| `-1` | `3` | hadronized hole |

Raw label `-2` is a hard-parton marker, not a final hadron. It is excluded
from all hadron vectors and all jet inputs. Its event count is retained in
`nHardMarkers`. Event-level counts separate normal, positive wake, negative
wake, label-2 negative thermal hadrons, and label-3 hadronized holes.
`signedPx`, `signedPy`, `signedPz`, `signedE`, and `signedSumPt` use a minus
sign for raw labels 2 and 3.

## Jet reconstruction

Jets are reconstructed independently in the two event variants with FastJet
anti-kT, E-scheme recombination, for `R=0.1`, `R=0.2`, `R=0.4`, and `R=0.8`.

- There is no constituent pT or eta cut. Zero-pT records and raw label `-2`
  markers are excluded.
- Normal and positive-wake hadrons enter clustering with their physical
  four-momenta.
- Each negative wake/hole is inserted into anti-kT clustering as an
  infinitesimal particle at the same rapidity and azimuth. Its physical
  four-vector is then subtracted from the jet that captures that ghost. This
  is the jet-level 4MomSub convention.
- Jets are stored when their pre-subtraction `RawPt` is at least 1 GeV and
  `abs(RawEta) < 5`. These are loose storage cuts; physics selections should be
  applied later.
- Stored jets remain ordered by pre-subtraction `RawPt`.

The primary `jet4Pt`, `jet4Eta`, `jet4Phi`, and `jet4mass` branches are the
four-momentum-subtracted values. A negative `jet4mass` means that subtraction
made the corrected four-vector spacelike; it is deliberately preserved rather
than clipped. `jet4RawPt`, `jet4RawEta`, `jet4RawPhi`, and `jet4RawMass`
retain the positive-constituent jet before subtraction. The same convention
applies to every `jet1*`, `jet2*`, and `jet8*` branch.

The nonlinear substructure observables are calculated from normal plus
positive-wake constituents. Negative contributions cannot be inserted into a
nonlinear clustering tree as ordinary particles. This scope is explicit in
the metadata and must not be described as constituent-level wake subtraction.
Per-jet `NegativeWakePt`, `NNegativeWake`, `NNegativeThermal`,
`NHadronizedHoles`, `PositiveWakePt`, and `WakeFraction` branches make the
effect auditable. `SignedG` is also provided as a signed linear diagnostic;
the requested `G` is the positive-constituent girth.

Schema v6 adds `TotalMult`, the signed count
`NNormal + NPositiveWake - NNegativeWake`. Here `NNegativeWake` includes both
raw-label-2 negative thermal hadrons and raw-label-3 hadronized holes. This is
a signed bookkeeping observable, not the multiplicity of a physical
negative-subtracted constituent list.

## Jet branches

Replace `X` by `1`, `2`, `4`, or `8` for R=0.1, R=0.2, R=0.4, or R=0.8:

| Branch | Definition |
| --- | --- |
| `jetXEta`, `jetXPhi`, `jetXPt`, `jetXmass` | 4MomSub-corrected jet kinematics |
| `jetXZg`, `jetXRg` | first Soft Drop passing split |
| `jetXMult` | number of normal plus positive-wake constituents |
| `jetXTotalMult` | signed count `NNormal + NPositiveWake - NNegativeWake` |
| `jetXPtD` | `sqrt(sum_i pT_i^2) / sum_i pT_i` |
| `jetXEffectiveMultiplicity` | `1/jetXPtD^2` for normal plus positive-wake constituents |
| `jetXLeadingFraction` | leading positive-constituent pT divided by positive scalar pT |
| `jetXNormalPtD` | momentum dispersion using raw-label-0 hadrons only |
| `jetXNormalEffectiveMultiplicity` | `1/jetXNormalPtD^2`, excluding all wake hadrons |
| `jetXLeadingNormalFraction` | leading raw-label-0 hadron pT divided by normal scalar pT |
| `jetXG` | `sum_i pT_i DeltaR(i,jet) / sum_i pT_i` |
| `jetXMaxKt` | largest `min(pT1,pT2) DeltaR12` in the full C/A tree |
| `jetXHardPartonId`, `jetXHardPartonPt`, `jetXHardPartonDR` | one-to-one matched outgoing hard-parton marker truth tag |

Soft Drop reclusters the positive constituents with Cambridge/Aachen and
follows the harder-pT branch. The parameters are `z_cut=0.1`, `beta=0`, and
`R0=R`. The stored quantities are

```text
zg = min(pT1,pT2)/(pT1+pT2),  Rg = DeltaR12.
```

`SoftDropValid=0` and floating-point `NaN` values mark a jet with no passing
split. Additional branches retain the groomed pT and mass, the number of
passing harder-branch splittings, and the z and angle of the maximum-kT split.
`MaxKtPrimary` gives the maximum along only the harder-branch declustering;
the requested `MaxKt` searches the full binary tree.

Jets in the paired variants are matched one-to-one by increasing raw-axis
distance with `DeltaR < R/2`. `PairMatchIndex`, `PairMatchDR`, and
`PairMatchOtherPt` permit direct no/with-prehydro comparisons.
`PairMatchOtherHardPartonId` records the marker tag on the matched jet so that
flavor-tagged analyses can require the same outgoing marker in both variants.
`PairMatchIndex=-1` means no match.

Outgoing hard-parton markers are matched independently to each variant by
increasing raw-axis distance with `DeltaR < R`. Each marker and jet can be used
only once. `HardPartonId=0` and non-finite pT/distance mean that no marker was
matched. Absolute PDG IDs 1--6 identify a quark tag and ID 21 identifies a
gluon tag. This is generator truth, not an experimentally accessible flavor
definition.

Schema v3 adds the normal-only fragmentation and hard-parton-match branches.
Schema v4 adds the complete `jet1*` branch family and corresponding event,
pair, and conversion totals for R=0.1. The weighted jet spectra, nonlinear
substructure overlays, migration-safe paired Soft Drop and momentum-dispersion
audit, pT-slice closure, and inclusive jet-RAA workflow all include R=0.1,
R=0.2, R=0.4, and R=0.8. The effective-charge response remains scoped to
R=0.2, R=0.4, and R=0.8.

## Formation-time estimator (schema v5, retained in schema v6)

Schema v5 stores the Cambridge--Aachen declustering tree needed for the
final-state formation-time estimator for R=0.4 and R=0.8 jets with
4MomSub-corrected `Pt > 30 GeV`. The calculation is performed in the C++ ROOT
writer from the original double-precision hadron four-vectors and FastJet
history. Reconstructing the tree later from the rounded float hadron branches
is not equivalent and is not supported for this observable.

For every valid internal node,

```text
tau_f [fm/c] = 0.19732698
               / (2 E_parent z1 z2 [1-cos(theta12)])
zi           = Ei / E_parent
z            = min(z1,z2)
kT           = min(pT1,pT2) DeltaR12.
```

`theta12` is the exact three-dimensional opening angle calculated from the
two child momentum vectors, and `E_parent` is the C/A parent energy. The
stored small-angle audit replaces `1-cos(theta12)` by `DeltaR12^2/2`; it is a
validation quantity, not the nominal estimator.

The C/A tree contains normal and positive-wake hadrons at physical
four-momentum. Raw-label-2 negative wake particles and raw-label-3 hadronized
holes are ghost-associated during anti-kT reconstruction and subtracted from
the jet four-vector. They therefore affect the corrected-pT selection but are
not inserted into the nonlinear C/A constituent tree. This is deliberate:
4MomSub does not define a unique negative-subtracted nonlinear tree.

For `jet4*` and `jet8*`, schema v5 adds the following branches, retained
unchanged in schema v6:

| Branch suffix | Definition |
| --- | --- |
| `FormationTauF` | flat vector of exact `tau_f` values for all valid full-tree splits |
| `FormationTauFSmallAngle` | corresponding DeltaR small-angle audit values |
| `FormationZ` | `min(E1,E2)/E_parent` |
| `FormationTheta` | exact three-dimensional opening angle |
| `FormationDeltaR` | rapidity-azimuth distance between the children |
| `FormationKt` | `min(pT1,pT2) DeltaR12` |
| `FormationParentE` | parent energy in GeV |
| `FormationOffset` | cumulative split offset for each stored jet; length is `nJet+1` |
| `FormationInvalidSplits` | per-jet count of nodes rejected by finite/positive checks |
| `FormationHardestValid` | one when a valid full-tree split exists |
| `FormationHardest*` | values at the global maximum-`kT` full-tree split |

The seven all-split branches are flat event vectors to remain compatible with
ROOT collection dictionaries. Splits for jet `j` occupy
`FormationOffset[j]:FormationOffset[j+1]`. For every retained jet, the writer
and analyzer require `valid splits + invalid splits = Mult - 1`. Jets outside
the declared radius/pT scope have empty split ranges and invalid count zero.

`analysis/plot_oo_jet_formation_time.py` produces event-weighted
`dN/dlog10(tau_f)` spectra, cross-section tables, exact/small-angle audits,
and `tau_f` versus `z`, `DeltaR`, and `kT` distributions. It reports every
valid full-tree declustering and, separately, one global hardest-`kT` split
per jet. Corrected-pT intervals are `(30,50]`, `(50,80]`, and `(80,infinity)`.
The biased-PYTHIA event weight is applied once. Ratios and integrated
diagnostics use paired delete-one-run jackknives; plotted one-dimensional
ratios require at least 20 no-prehydro entries in a bin.

For the three-way V3 comparison, run the pair analyzer once on each aligned
ROOT file and merge with `analysis/plot_oo_v3_jet_formation_time.py`. The
merger refuses differing no-prehydro spectra, summaries, or two-dimensional
tables. Cambridge--Aachen declustering here is a formation-time estimator
reconstructed from final jet constituents. It does **not** reproduce the
generator-level parton-shower history or identify the actual splitting time
of a shower parton.

The preferred final-state effective-charge proxy is

```text
N_eff_normal = (sum_normal pT)^2 / sum_normal(pT^2)
             = 1 / NormalPtD^2.
```

It is one for a single normal constituent and grows when momentum is shared
among many normal constituents. Excluding raw labels 1--3 prevents the wake
from manufacturing apparent shower charges. It is still a final-state
hadron-level proxy: it does not equal the number of active or medium-resolved
shower partons at the hydro start.

## Running

Use a fixed local snapshot for every source. Do not point a conversion at a
directory that a production monitor is actively updating.

```bash
python3 production/oo_planb_2509/analysis/convert_oo_paired_to_root.py \
  --source v2=/path/to/frozen/v2 \
  --aa-task-manifest /path/to/aa_task_manifest.tsv \
  --output /path/to/oo_paired.root \
  --build-dir /raid5/data/yjlee/hybrid_dev/test/tmp_oo_root_build
```

The ROOT file embeds the source mapping, cuts, schema version, and conversion
totals. The JSON summary adds the ROOT and FastJet versions, source-repository
commit, command, output size, and SHA256. The converter validates all required
trees, entry counts, and branches with uproot before atomically publishing the
output. When supplied, the manifest path, row count, and SHA256 are also
recorded in the JSON summary; accepted rows in the audit TSV retain the hydro
event ID, Ncoll, and payload SHA256.

For the 50,000-event v2 campaign,
`analysis/run_oo_v2_final_analysis.sh` guards this conversion on the strict
completion marker, runs the complete-prefix hadron RAA analysis, creates the
manifest-validated ROOT file, produces the inclusive and four-slice jet
comparisons, and requires exact 50,000-pair acceptance and slice closure within
`1e-12` mb before writing its own completion marker.

Run the parser tests with:

```bash
python3 production/oo_planb_2509/analysis/test_convert_oo_paired_to_root.py -v
python3 production/oo_planb_2509/analysis/test_plot_oo_jet_formation_time.py -v
```

For the weighted no/with-prehydro jet-variable comparison, use
`analysis/plot_oo_jet_variables.py`. It applies the one-event-run
`PythiaParallel` normalization once, selects on the corrected jet pT, and uses
a paired delete-one-run jackknife for distributions and ratios. Absolute
cross sections remain in the TSV audit, while plotted shapes use
`(1/sigmaJetSelected) dSigma/dx` separately for each variant. Shape ratios
include that selected-jet normalization in every jackknife replica. Its weighting
regression test is `analysis/test_plot_oo_jet_variables.py`. An optional
`--pt-max` is inclusive while `--pt-min` remains exclusive. Consequently the
four nonoverlapping campaign intervals are generated with `(20,30]`,
`(30,50]`, `(50,80]`, and `(80,infinity)`. The companion
`analysis/summarize_oo_jet_pt_slices.py` verifies their normalization metadata
and exact cross-section closure to the inclusive `pT > 20` result before
writing the combined table and figure.

The `Zg` and `Rg` histograms include all selected jets. Bin zero is an
explicit `softdrop_failed` sentinel bin for `SoftDropValid=0`; its width is
one physical bin and its integral is therefore the failed-jet cross section.
For `Zg` the sentinel interval is `[0.075,0.1)`. For `Rg` it is the
equal-width interval immediately below zero, with a radius-dependent width.
All remaining bins contain `SoftDropValid=1` jets. The histogram TSV records
the bin role in `bin_kind`, while the summary and JSON metadata retain raw and
weighted pass/fail counts and fractions. Reported `Zg` and `Rg` moments remain
conditioned on `SoftDropValid=1` and never average in the sentinel value.

For the matched effective-charge sensitivity study, use
`analysis/plot_oo_jet_charge_response.py`. It selects and bins on the
no-prehydro jet, requires its one-to-one Plan-B match, and measures
`1 - pT(Plan B)/pT(no prehydro)`. Positive values mean additional loss from
turning on prehydro. Results are provided versus wake-excluded effective
multiplicity, leading-normal fraction, `NSD+1`, and maximum-kT splitting for
all, quark-tagged, and gluon-tagged jets. Quark/gluon subsets require identical
hard-marker PDG IDs in both variants. The nominal pair-axis requirement is
`DeltaR<R/2`; `--max-pair-match-dr-fraction` supports stricter robustness
checks. Conditional means use the biased PYTHIA event weight once, and
uncertainties use a paired delete-one-run jackknife.

For a migration-safe comparison of Soft Drop and momentum dispersion, use
`analysis/plot_oo_jet_paired_substructure.py`. The corrected pT interval and
eta acceptance are imposed only on the no-prehydro jet. Its one-to-one
`PairMatchIndex` partner is retained even if Plan B moves it outside that pT
interval. The output stores the four no-prehydro-to-Plan-B Soft Drop states,
paired `PtD` and `NormalPtD` shifts (including common-sample ratios of means),
positive and signed multiplicities, both-pass `Zg` and `Rg`, girth, maximum
`kT`, and R=0.4/0.8 all-split and hardest-split formation-time summaries.
Formation-time rows report the arithmetic per-jet mean, mean `log10(tau_f)`,
and hardest-`kT` `log10(tau_f)` response. The output also retains incremental
fractional loss for all,
quark-tagged, gluon-tagged, Soft-Drop, and wake-excluded effective-multiplicity
categories. Event-level jackknife errors are the nominal statistical errors.
The delete-one-`hydroIndex` block error, effective number of weighted event
contributors, and largest single-event weight fraction provide correlated
hydro and weight-tail checks.

Because `NormalEffectiveMultiplicity=1/NormalPtD^2`, conditioning on its
single-like or many-like ranges and then inspecting a `NormalPtD` shift has a
built-in boundary correlation. Those rows diagnose migration; the matched-pT
loss versus the fixed no-prehydro class is the corresponding physics-facing
test.

The negative-particle treatment follows the jet-level ghost-association and
four-vector-subtraction construction described for 4MomSub in
[arXiv:1612.05116](https://arxiv.org/abs/1612.05116). Its use for the negative
Hybrid wake is discussed in
[JHEP 01 (2025) 164](https://arxiv.org/abs/2409.12238). The nonlinear-shape
caveat above remains important: jet-level 4MomSub defines corrected jet
four-momentum, not a unique corrected constituent tree.

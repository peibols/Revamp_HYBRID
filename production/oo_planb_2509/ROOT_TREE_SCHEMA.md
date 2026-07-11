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
in the audit TSV.

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
anti-kT, E-scheme recombination, for `R=0.4` and `R=0.8`.

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
applies to every `jet8*` branch.

The nonlinear substructure observables are calculated from normal plus
positive-wake constituents. Negative contributions cannot be inserted into a
nonlinear clustering tree as ordinary particles. This scope is explicit in
the metadata and must not be described as constituent-level wake subtraction.
Per-jet `NegativeWakePt`, `NNegativeWake`, `NNegativeThermal`,
`NHadronizedHoles`, `PositiveWakePt`, and `WakeFraction` branches make the
effect auditable. `SignedG` is also provided as a signed linear diagnostic;
the requested `G` is the positive-constituent girth.

## Jet branches

Replace `X` by `4` or `8` for R=0.4 or R=0.8:

| Branch | Definition |
| --- | --- |
| `jetXEta`, `jetXPhi`, `jetXPt`, `jetXmass` | 4MomSub-corrected jet kinematics |
| `jetXZg`, `jetXRg` | first Soft Drop passing split |
| `jetXMult` | number of normal plus positive-wake constituents |
| `jetXPtD` | `sqrt(sum_i pT_i^2) / sum_i pT_i` |
| `jetXG` | `sum_i pT_i DeltaR(i,jet) / sum_i pT_i` |
| `jetXMaxKt` | largest `min(pT1,pT2) DeltaR12` in the full C/A tree |

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
`PairMatchOtherPt` permit direct no/with-prehydro comparisons. `-1` means no
match.

## Running

Use a fixed local snapshot for every source. Do not point a conversion at a
directory that a production monitor is actively updating.

```bash
python3 production/oo_planb_2509/analysis/convert_oo_paired_to_root.py \
  --source original=/path/to/frozen/original \
  --source continuation=/path/to/frozen/continuation \
  --output /path/to/oo_paired.root \
  --build-dir /raid5/data/yjlee/hybrid_dev/test/tmp_oo_root_build
```

The ROOT file embeds the source mapping, cuts, schema version, and conversion
totals. The JSON summary adds the ROOT and FastJet versions, source-repository
commit, command, output size, and SHA256. The converter validates all required
trees, entry counts, and branches with uproot before atomically publishing the
output.

Run the parser tests with:

```bash
python3 production/oo_planb_2509/analysis/test_convert_oo_paired_to_root.py -v
```

For the weighted no/with-prehydro jet-variable comparison, use
`analysis/plot_oo_jet_variables.py`. It applies the one-event-run
`PythiaParallel` normalization once, selects on the corrected jet pT, and uses
a paired delete-one-run jackknife for distributions and ratios. Its weighting
regression test is `analysis/test_plot_oo_jet_variables.py`. An optional
`--pt-max` is inclusive while `--pt-min` remains exclusive. Consequently the
three nonoverlapping campaign intervals are generated with `(30,50]`,
`(50,80]`, and `(80,infinity)`. The companion
`analysis/summarize_oo_jet_pt_slices.py` verifies their normalization metadata
and exact cross-section closure to the inclusive `pT > 30` result before
writing the combined table and figure.

The negative-particle treatment follows the jet-level ghost-association and
four-vector-subtraction construction described for 4MomSub in
[arXiv:1612.05116](https://arxiv.org/abs/1612.05116). Its use for the negative
Hybrid wake is discussed in
[JHEP 01 (2025) 164](https://arxiv.org/abs/2409.12238). The nonlinear-shape
caveat above remains important: jet-level 4MomSub defines corrected jet
four-momentum, not a unique corrected constituent tree.

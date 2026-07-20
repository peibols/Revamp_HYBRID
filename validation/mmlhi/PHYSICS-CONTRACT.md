# MMLHI Heavy-Transport Physics Contract

This document separates the backward-compatible runtime defaults from the
configuration used for controlled heavy-flavor validation. It is a diagnostic
transport contract, not a complete matching prescription or heavy-flavor tune.

## No-Overlap Validation Configuration

```text
heavy_quark_eloss_mode = 2
heavy_quark_lambda = 1.961
heavy_quark_charm_mass = 1.25
heavy_quark_bottom_mass = 4.2
heavy_quark_add_generic_broadening_with_diffusion = false
heavy_quark_enable_hard_moliere = false
```

Mode 2 compares heavy drag with the light-HYBRID strong-coupling loss in the
local fluid frame. It applies the smaller allowed energy loss while preserving
the heavy-particle mass shell. As in the legacy implementation, heavy diffusion
is added only when the heavy-drag branch wins.

The configured masses are lower invariant-mass floors. They do not overwrite a
larger live shower mass. The pinned PYTHIA 8.315 particle data uses
`m0(c)=1.50 GeV` and `m0(b)=4.80 GeV`, so direct heavy quarks normally retain
those larger masses. The per-step increment is
`dt_star = gamma_flow*(1-v_flow dot v_parton)*dt`, the elapsed time in the local
fluid frame. The C++ field retains the historical `fluid_path_length_fm` name;
calling it a spatial distance is only accurate in the ultrarelativistic limit.

The two Boolean settings default to `true` to preserve results from the first
MMLHI implementation when they are omitted. They are set to `false` in the
validation template for the reasons below.

## Soft Momentum-Transfer Matching

`heavy_quark_add_generic_broadening_with_diffusion = false` suppresses the
generic light-parton Gaussian `kappa` kick for heavy quarks in modes 2 and 3.
Light partons retain the existing generic broadening. Mode 1 has no heavy
diffusion and therefore retains generic broadening independently of this
switch.

There is an important limitation. In Mode 2, heavy diffusion is coupled to the
heavy-drag branch. If the light-HYBRID baseline wins, the current no-overlap
setting applies neither generic broadening nor heavy diffusion on that step.
The alternative compatibility setting applies generic broadening on every
step and heavy diffusion on drag-winning steps, so those drag steps contain
two additive soft descriptions. The code now exposes both choices, but neither
is claimed as a finished hard/soft matching prescription.

The no-overlap setting is useful for implementation comparisons because it
does not add two soft kicks. Production still needs an explicit decision:
apply heavy diffusion also on baseline-winning steps, retain generic broadening
only on those steps, or derive a matched subtraction. The Moliere
implementation itself uses a soft Gaussian sector below a hard-transfer
threshold and rare explicit hard scatterings above it; the threshold is
intended to control double counting (arXiv:2603.08776).

## Hard Moliere Scattering

`heavy_quark_enable_hard_moliere = false` disables explicit hard Moliere
sampling for charm while leaving light-quark and gluon Moliere scattering
unchanged. Bottom is always excluded because the existing hard-scattering
tables and flavor-matching path only cover flavors through charm.

The existing charm path is retained behind the switch for controlled software
tests, but it is not recommended for physics production. A 100-event forced
charm audit found seven accepted hard charm scatterings. All seven failed both
the outgoing heavy-mass-shell check and immediate
`p_projectile + p_thermal = p_projectile' + p_recoil` closure. The largest
residuals were 2.25 GeV^2 in mass squared and 0.447 GeV in the Euclidean
four-vector closure norm. These failures arise because the inherited table
sampler constructs massless projectile kinematics.

Projecting only the outgoing charm back onto a massive shell is not an
acceptable repair: it would not make the sampled rate, angular distribution,
or recoil kinematics into a consistent massive `2 -> 2` process. Hard charm
should remain off until the sampler and tables use massive incoming and
outgoing kinematics end to end.

## Bottom Policy

Bottom receives the continuous mode-1/2/3 heavy kernel. It does not receive
explicit hard Moliere scattering. Setting `do_elastic = true` still selects
the existing Moliere propagation and hadronization path for the rest of the
event, so parton-level and hadron-level comparisons must distinguish that
global selection from direct scattering of the bottom quark.

## Finite LRES Composition

The heavy kernel acts on the effective object selected by the existing LRES
timeline and Moliere Mode A-E algorithm. It does not change the LRES resolution
clock. A heavy daughter receives heavy transport only while that daughter is a
materialized active object; an unresolved light parent receives light/coherent
transport. This mode-dependent identity assignment is a physics choice and
must be included in any Mode A-E interpretation.

Mode D daughter probes do not commit momentum, recoil, hole, or heavy-step
diagnostics, but they do advance copied RNG state in a fixed daughter-1 then
daughter-2 sequence and return that state to the global transport generator.
Only Mode E uses deterministic branch-local probe streams and restores the
global elastic generator after the probe batch. “Side-effect-free” must
therefore not be used as an unqualified description of Mode D.

There are two additional Mode-E qualifications:

- The current live-parent opening preserves total daughter energy but not the
  full parent three-momentum. In the 20-event feature-on sample, the 23 opening
  operations have average/max spatial residuals of `0.189/1.856 GeV`; the
  maximum relative residual is `0.222`. Energy residuals are zero at printed
  precision. This is not four-momentum conservation.
- The opening helper creates massless daughters. A seed-`870001` event-display
  audit materializes an anti-charm daughter massless at `t=6.166 fm/c`; its
  first in-medium heavy update at `t=7.166 fm/c` imposes the `1.25 GeV` floor
  and changes lab energy from `48.768` to `51.351 GeV`. A mass-aware opening
  prescription is required before Mode E is used for heavy-flavor physics.

The feature-on 20-event C/D/E runs sampled no unresolved elastic candidate.
They establish continuous heavy execution through those paths, while the
resolving-scattering and response-accounting tests are currently heavy-mode-off.

## Hadronization Interface

The current `LundGenerator` keeps each transported three-momentum but rebuilds
its energy from PYTHIA's particle-data mass before appending the parton. It does
not pass the transport energy through unchanged. This re-on-shell operation is
not an integrated heavy-hadron formation model and can add another energy
discontinuity when a transported heavy daughter is on the configured floor
rather than PYTHIA's `m0`. No D/B-hadron claim should use this interface without
an explicit fragmentation/coalescence and four-momentum-closure validation.

## Remaining Physics Gates

The no-overlap configuration is suitable for implementation validation and
controlled parton-level comparisons. Production D- or B-hadron claims still
require:

- a chosen heavy-hadron formation contract, including fragmentation and any
  intended recombination/coalescence;
- sufficient-statistics charm and bottom spectra in the target medium;
- observable-level `R_AA` and `v2` validation;
- a massive hard-scattering implementation before enabling hard Moliere for
  charm;
- a mass-aware, four-momentum-defined Mode-E daughter materialization before
  combining recursive coherence with heavy transport;
- an explicit soft-matching choice for baseline-winning Mode-2 steps;
- a separate physics decision and tables before adding hard Moliere for
  bottom.

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
- an explicit soft-matching choice for baseline-winning Mode-2 steps;
- a separate physics decision and tables before adding hard Moliere for
  bottom.

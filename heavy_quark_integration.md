# Heavy-Quark Energy-Loss Integration

## Scope

This branch layers a heavy-quark transport kernel on
`main_moliere_lres_integration`. It integrates heavy-quark drag and diffusion
into the existing MMLI propagation paths. It does not yet port the legacy
heavy-flavor recombination, coalescence, or dedicated heavy-hadron analysis.

The implementation is default-off. With `heavy_quark_eloss_mode = 0`, it does
no heavy-specific arithmetic and consumes no additional random numbers.

## Runtime Interface

```text
heavy_quark_eloss_mode = 0
heavy_quark_lambda = 1.961
heavy_quark_charm_mass = 1.25
heavy_quark_bottom_mass = 4.2
heavy_quark_add_generic_broadening_with_diffusion = true
heavy_quark_enable_hard_moliere = true
```

`heavy_quark_lambda` must be supplied explicitly when a nonzero mode is
selected. The legacy code convention is

```text
kappa_HQ = pi * sqrt(heavy_quark_lambda).
```

For example, `kappa_HQ = 4.4` corresponds to
`heavy_quark_lambda = (4.4/pi)^2`, approximately `1.96`. The mass settings are
floors. If a shower parton already has a larger timelike invariant mass, the
kernel preserves that larger mass.

The two Boolean options default to `true` so configurations from the first
MMLHI implementation retain their behavior. The controlled no-overlap
validation setup sets both to `false`: it suppresses generic broadening in
modes 2 and 3, and hard charm Moliere remains disabled until the inherited
massless scattering sampler is replaced. This is not yet a complete soft
matching prescription; see the explicit
[`validation/mmlhi/PHYSICS-CONTRACT.md`](validation/mmlhi/PHYSICS-CONTRACT.md).

Modes 1 and 2 require the global light-parton setting `mode = 0`:

| Heavy mode | Behavior for charm and bottom |
|---|---|
| 0 | Disabled; exact parent-MMLI path |
| 1 | Choose the smaller loss between heavy drag and the light-HYBRID strong-coupling candidate |
| 2 | Mode 1, with heavy diffusion added whenever the heavy-drag branch is selected |
| 3 | Heavy diffusion only; the previously advertised but broken legacy mode is now explicit |

## Per-Step Algorithm

For a charm or bottom quark in a cell with `T >= Tc`:

1. The existing MMLI soft Gaussian broadening is applied if `kappa != 0` and
   the soft-matching switch allows it. It is always retained for light partons
   and mode 1; modes 2 and 3 can use heavy diffusion instead.
2. The four-momentum is boosted to the local fluid rest frame.
3. The effective heavy mass is the larger of the configured floor and the
   current timelike invariant mass.
4. For modes 1 and 2, form the drag coefficient
   `eta_D = (pi/2) sqrt(lambda) T^2/M` and apply the first-order spatial update
   `p_i -> p_i (1 - eta_D*dx/0.2)`. Energy is recomputed on shell.
5. Compare that drag energy loss to the light-HYBRID strong-coupling candidate
   in the same fluid frame. Use heavy drag when it loses less energy, when the
   light candidate would cross the mass shell, or when the light stopping
   distance is exhausted. Otherwise, apply the winning baseline energy loss in
   the same fluid frame while preserving the heavy mass shell.
6. In mode 2, add independent Gaussian kicks with per-component variance
   `pi*sqrt(lambda)*gamma*T^3*dx/0.2` after a selected drag update.
7. In mode 3, apply only the diffusion update.
8. Put the heavy quark on shell and boost back to the lab frame.

The factor `1/0.2` preserves the codebase's existing GeV/fm conversion. The
new on-shell reconstruction intentionally removes the order-step-squared
energy inconsistency in the legacy Euler update. Feature-on random sequences
therefore are not expected to be byte-identical to the old heavy executable.
Diffusion trajectories differ even with aligned generator seeds because the
legacy Box-Muller draw uses a sine phase while this kernel uses a cosine phase.
The distribution-level contract is the shared zero mean and variance, which is
tested directly.

## Composition With MMLI

- **Standard propagation:** the kernel runs inside `EnergyLoss::loss_rate`.
- **Finite LRES:** continuous loss acts on whichever effective object the
  existing LRES timeline makes active. A heavy daughter receives heavy
  transport only when the active object itself has charm or bottom identity.
- **Moliere modes A-E:** the same kernel runs inside the Moliere segment
  stepper. Daughter-level segments receive daughter-level heavy transport.
  Mode D/E speculative candidate probes include the heavy kinematics but do
  not increment committed-step diagnostics.
- **Hard Moliere plus diffusion:** eligible hard scattering occurs before the
  continuous heavy update. Charm hard scattering is controlled separately by
  `heavy_quark_enable_hard_moliere`; bottom is always excluded. Generic soft
  broadening is additive only when its matching switch is enabled.
- **Soft-matching caveat:** Mode-2 heavy diffusion runs only when the drag
  branch wins. With generic broadening disabled, a baseline-winning step has
  no stochastic soft kick. With the compatibility setting enabled, a
  drag-winning step receives both generic and heavy diffusion. Production
  needs an explicit prescription between these two diagnostic limits.
- **Recoilers:** a charm recoiler propagated by the existing rescattering loop
  also receives the heavy kernel.
- **Bottom:** bottom receives continuous heavy drag/diffusion. The current
  Moliere hard-scattering tables and matching logic are charm-specific, so
  bottom remains excluded from hard Moliere sampling.
- **Sources/wakes:** no separate heavy response model is added. Builds with
  the existing source bookkeeping see the net step-level momentum change.

## Diagnostics

When a heavy mode is enabled, the `EnergyLoss` destructor reports:

```text
n_heavy_steps
n_baseline_steps
n_drag_steps
n_diffusion_steps
n_diffusion_only_steps
n_invalid_steps
sum_energy_change
n_hard_scattered_heavy
n_hard_heavy_mass_shell_failures
n_hard_heavy_momentum_closure_failures
max_hard_heavy_mass_shell_residual
max_hard_heavy_momentum_closure_residual
```

`sum_energy_change` is the signed lab-frame energy change from committed heavy
updates. It may be negative for diffusion-only trajectories because stochastic
kicks can add energy. The hard-scattering counters inspect the immediate
accepted `2 -> 2` record before later soft broadening or continuous loss.

## Validation Performed

Validation used PYTHIA 8.315. The reproducible runner and exact results are in
[`validation/mmlhi`](validation/mmlhi/README.md) and
[`RESULTS-20260720.md`](validation/mmlhi/RESULTS-20260720.md).

- Standalone deterministic tests cover disabled-mode RNG closure, light-parton
  bypass, analytic drag, crossover selection, diffusion determinism, explicit
  mode 3, moving fluid, charm/bottom mass shells, and invalid configuration.
- Parent MMLI and MMLHI mode 0 are byte-identical for one-event standard,
  finite-LRES, Moliere, and finite-LRES-plus-Moliere runs.
- A forced low-`pThat` charm event exercises modes 1, 2, and 3; mode 2 reruns
  byte-identically.
- Forced charm finite-LRES, Moliere, and combined finite-LRES-plus-Moliere
  smokes complete with nonzero committed heavy-step counters.
- Heavy-mode-off Modes A-E close byte-for-byte to the parent MMLI branch.
  Targeted C, D, and E samples exercise resolving scatterings; candidate
  accounting closes, Mode E has zero frontier-order mismatches, and each
  accepted unresolved scattering produces exactly one recoil and one hole.
- A controlled same-tree comparison to the clean legacy heavy branch gives
  small parton-level differences at `|eta| < 1`: 3.64% RMS in charm `pT` and
  1.58% RMS in bottom `pT` across 11 matched final heavy quarks of each flavor.
- A 100-event forced-charm Moliere audit accepts seven hard charm scatterings.
  All seven expose the known massless-table gap: they fail both the heavy mass
  shell and immediate projectile-plus-medium four-momentum closure checks.
- A forced bottom event completes in standard and Moliere configurations. The
  parton outputs agree when no hard scattering is sampled, as expected from
  the current bottom exclusion. Hadron files differ because requesting
  Moliere selects the existing Moliere hadronization path.

These are implementation and smoke validations, not sufficient-statistics
physics validation. The no-overlap configuration prevents generic soft
`kappa` broadening from being added to heavy diffusion and disables hard charm
Moliere, but leaves baseline-winning Mode-2 steps without a stochastic soft
kick. A production claim still requires a resolved soft-matching prescription,
a massive hard-scattering implementation if hard charm is desired, plus the
intended D/B-hadron formation contract and observable-level `R_AA` and `v2`
validation.

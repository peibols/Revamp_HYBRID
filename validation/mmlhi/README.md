# MMLHI Validation Matrix

This directory contains the reproducible software-validation matrix for the
heavy-quark transport layer on top of MMLI. MMLHI inherits the complete MMLI
light-parton transport path; its intended specialization is that charm and
bottom remain massive transported objects instead of taking MMLI's early
free-stream return. Exact parent/child closure therefore uses a light-only
PYTHIA card. Heavy-containing runs must differ and are checked against the
massive-transport contract instead. The matrix uses the pinned PYTHIA 8.315
installation by default and writes all generated inputs, binaries, logs, and
event output under the workspace `test/` directory.
The callback-contract check uses the Python `uproot` package to inspect the
generated ROOT event-display tree.

Run from the MMLHI worktree:

```bash
validation/mmlhi/run_validation.sh
```

The latest validated matrix and its physics interpretation are recorded in
[`RESULTS-20260720.md`](RESULTS-20260720.md). The controlled runtime choices,
their matching caveats, and remaining production gates are separated in
[`PHYSICS-CONTRACT.md`](PHYSICS-CONTRACT.md).

The defaults assume the validated local MMLI parent, Moliere tables, PbPb smoke
hydro inputs, and the recorded high-pT Mode-E seed fixture in `clean_port_run`.
They can be overridden with `PARENT_REPO`, `PYTHIA_HOME`, `MOLIERE_TABLES`,
`MMLHI_REFERENCE_RUN`, `MMLHI_MODEE_REFERENCE_RUN`, and
`MMLHI_VALIDATION_OUT`.

The runner checks:

- pinned parent and child builds;
- the Mode-E coherent-source eligibility policy for quarks, gluons, photons,
  leptons, representative hadrons, and PDG ID zero;
- unit tests plus AddressSanitizer and UndefinedBehaviorSanitizer;
- byte-identical heavy-mode-off closure in standard, LRES, Moliere, and
  combined Mode A-E propagation;
- byte-identical standard Moliere output with a null scattering callback and
  an `Apply` observer callback, with ROOT records proving that hard candidates
  and one recoil/hole pair per candidate were observed;
- targeted Mode-C, Mode-D, and recursive Mode-E resolving-scattering paths;
- seed-fixed Mode-E failed-daughter veto, independent coherent-source accept,
  and resolving parent-candidate veto paths in the inherited MMLI core (the
  historical fixture contains charm and is not an exact heavy-branch oracle);
- forced charm modes 1, 2, and 3;
- forced charm with LRES, Moliere, and combined Modes A-E;
- bottom continuous transport with Moliere requested;
- deterministic mode-2 and Mode-E reruns;
- zero hard-heavy scatters in the production-validation cards; the inherited
  massless hard-heavy sampler remains outside the production contract;
- missing-`heavy_quark_lambda` rejection;
- output completeness, invalid-step counters, dynamic candidate accounting,
  one recoil/hole pair per accepted unresolved scattering, finite output
  values, Mode-E frontier-permutation diagnostics, and explicit counts of
  PYTHIA hadronization retries or give-ups;
- Mode-E opening energy, absolute spatial-momentum, and relative
  spatial-momentum residuals. These are reported rather than treated as closed:
  the current live-parent opening is energy-preserving, not a four-momentum
  conserving daughter materialization.

The current August 13 matrix is written to
`test/mmlhi_validation_20260813`. Its complete counters and hashes are in
`summary.tsv`.

The matrix is an implementation validation. It is not a heavy-flavor physics
tune and does not establish D- or B-hadron observables because the legacy
heavy hadronization/coalescence contract has not been integrated.

For controlled same-tree comparisons against a legacy heavy output, use:

```bash
validation/mmlhi/compare_legacy_partons.py \
  LEGACY_RESULT_PARTONS MMLHI_RESULT_PARTONS --max-abs-eta 1
```

To inspect heavy-parton mass-shell changes in a ROOT event display, use:

```bash
validation/mmlhi/audit_event_display_mass_shell.py eventDisplay.root
```

The July 20 seed-870001 Mode-E audit output is stored in
[`results/20260720-modee-heavy-mass-shell-audit.tsv`](results/20260720-modee-heavy-mass-shell-audit.tsv).
It captures the currently known massless-daughter-to-heavy-floor transition.

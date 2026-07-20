# MMLHI Validation Matrix

This directory contains the reproducible software-validation matrix for the
heavy-quark transport layer on top of MMLI. It uses the pinned PYTHIA 8.315
installation by default and writes all generated inputs, binaries, logs, and
event output under the workspace `test/` directory.

Run from the MMLHI worktree:

```bash
validation/mmlhi/run_validation.sh
```

The latest validated matrix and its physics interpretation are recorded in
[`RESULTS-20260719.md`](RESULTS-20260719.md).

The defaults assume the validated local MMLI parent, Moliere tables, and PbPb
smoke hydro inputs. They can be overridden with `PARENT_REPO`, `PYTHIA_HOME`,
`MOLIERE_TABLES`, `MMLHI_REFERENCE_RUN`, and `MMLHI_VALIDATION_OUT`.

The runner checks:

- pinned parent and child builds;
- unit tests plus AddressSanitizer and UndefinedBehaviorSanitizer;
- byte-identical heavy-mode-off closure in standard, LRES, Moliere, and
  combined Mode A-E propagation;
- targeted Mode-C, Mode-D, and recursive Mode-E resolving-scattering paths;
- forced charm modes 1, 2, and 3;
- forced charm with LRES, Moliere, and combined Modes A-E;
- a 100-event charm hard-scattering coverage sample;
- bottom continuous transport with Moliere requested;
- deterministic mode-2 and Mode-E reruns;
- missing-`heavy_quark_lambda` rejection;
- output completeness, invalid-step counters, dynamic candidate accounting,
  one recoil/hole pair per accepted unresolved scattering, finite output
  values, and Mode-E frontier-permutation diagnostics.

The matrix is an implementation validation. It is not a heavy-flavor physics
tune and does not establish D- or B-hadron observables because the legacy
heavy hadronization/coalescence contract has not been integrated.

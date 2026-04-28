# lres Random-Seed Validation Notes

## PYTHIA Version Rule

For `lres_code_randomseedFix` validation builds and reproducibility checks, use the same validated PYTHIA install used for the `main`, `main_moliere_integration`, and `moliere_code_randomseedFix` campaigns:

- include: `/data/yjlee/pythia/pythia8/pythia8315/include`
- lib: `/data/yjlee/pythia/pythia8/pythia8315/lib`
- share: `/data/yjlee/pythia/pythia8/pythia8315/share/Pythia8`

Do not use arbitrary system PYTHIA headers or libraries for validation. In particular, do not silently fall back to `/usr/include`, `/lib`, or other local installations.

## Validation Requirement

Before accepting an `lres` validation result:

1. build the executable against the pinned PYTHIA 8.315 install above
2. confirm the resulting binary links to `/data/yjlee/pythia/pythia8/pythia8315/lib/libpythia8.so`
3. only then run the reproducibility or usage-mode validation

## Scope

This rule applies to:

- `pythielsen`
- `photelsen`
- `plasma`
- `plasma_debug`
- `plasma_lund`
- `wake`
- `hadronizer`

The goal is to avoid mixing validated physics results with random build environments.

## Trigger Workflow

`lres` trigger validation is a staged workflow, not a `main_trigger` executable. The branch trigger runner in [`test/run_it_trigger.sh`](/raid5/data/yjlee/hybrid_dev/Revamp_HYBRID/test/run_it_trigger.sh) should execute:

1. `photelsen`
2. `plasma_lund`
3. `wake`
4. `hadronizer`

The branch test directory must also contain:

- `setup_pythia.cmnd`
- `hydroinfoPlaintxtHuichaoFormat.dat`
- `TAb2LL.dat`

## No-Quench Reference Choice

For `lres` to produce a usable unquenched reference for comparison with `main_moliere_integration`, zero couplings are now treated as an explicit no-quench reference mode inside:

- [`plasma.cc`](/raid5/data/yjlee/hybrid_dev/Revamp_HYBRID/plasma.cc)
- [`plasma_lund.cc`](/raid5/data/yjlee/hybrid_dev/Revamp_HYBRID/plasma_lund.cc)
- [`plasma_debug.cc`](/raid5/data/yjlee/hybrid_dev/Revamp_HYBRID/plasma_debug.cc)

This choice was needed because the old staged wrapper did not have a true vacuum bypass:

- it still sampled geometry
- it still rewrote the tree through the staged medium wrapper
- and with `mode = 0` it could zero parton momenta in hot cells because the external `lres` convention (`0 = strong`) did not match the internal plasma-kernel convention (`3 = strong`)

The implemented rule is:

1. keep the external `lres` interface as `mode = 0` for strong coupling
2. normalize that to the internal legacy kernel value `qmethod = 3`
3. if `alpha = 0` and `kappa = 0`, run an explicit no-quench reference path:
   - preserve outgoing parton momenta
   - force `X = 0`, `Y = 0`

This fixes the invalid zero-momentum reference behavior. It does not by itself make `lres` exactly identical to `main_moliere_integration`; the staged wrapper still has different hard-header / remnant ordering and a different hadronization contract.

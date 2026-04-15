# randomseedFix

This branch is reserved for the random-seed and reproducibility cleanup in the HYBRID revamp.

## Goal

Make the simulation reproducible in a controlled and documented way by removing mixed RNG usage, making all random streams explicit, and verifying that the full quenched path reproduces byte-for-byte.

## Implemented Fixes

The seed cleanup is now implemented in the config-driven `main` path.

- Added an explicit `seed_base` config option. If it is not provided, the code falls back to `njob` for backward compatibility.
- Derived three explicit streams from that base:
  - `shower_seed = seed_base + 33`
  - `hybrid_seed = seed_base + 1346`
  - `lund_seed = seed_base + 2337`
- Threaded the explicit shower seed into `TreeGenerator`.
- Threaded the explicit Lund seed into `LundGenerator`.
- Printed the effective seeds at startup so the run provenance is visible in the log.
- Removed the stray `rand()` usage from the `ebe_hydro = 1` IPSAT path.
- Fixed the IPSAT sampling call so the sampler receives the full `Ncollsize_` range rather than a pre-drawn random integer.

## Remaining Seeding Bug That Was Found

After the first seed unification pass, the full

- `do_quench = true`
- `do_wake = true`
- `ebe_hydro = 0`

path still failed a same-seed reproducibility test.

The hard event weights and cross sections matched between reruns, but the production points and later wake/hadron outputs diverged starting from event 1. This showed that the explicit seed wiring was not the last source of nondeterminism.

The actual remaining bug was in `WakeGenerator`. The wake Metropolis loop used wall-clock time through `clock()` to decide when to restart or force a residual wake particle. That means two runs with the same seed could take different branches under different machine load, which changed the accepted wake configuration and then changed the later consumption of the shared HYBRID RNG stream.

In other words, the code was seeded, but it was not deterministic because one algorithm still depended on runtime timing.

## Final Wake Fix

The wake generator has now been made deterministic by replacing the wall-clock bailout with a proposal-count budget.

- Removed the `clock()` / `CLOCKS_PER_SEC` logic from `WakeGenerator`.
- Replaced it with a deterministic proposal budget:
  - base budget `200000`
  - restart increment `20000`
  - hard restart cap `100`
- Added a shared helper to add the residual wake particle when the remaining four-momentum is close enough to a physical hadron mass shell.

This preserves the original intent of the safeguard against pathological convergence, but makes the control flow depend only on the seeded RNG stream and the event state, not on machine timing.

## Verified Behavior

The fix was tested on the full averaged-hydro quenched path with wake generation enabled.

Test configuration:

- `Nev = 10`
- `do_quench = true`
- `do_wake = true`
- `ebe_hydro = 0`
- same hydro and Glauber inputs
- same executable and PYTHIA build

Results:

- Two clean runs with `seed_base = 0` produced byte-identical outputs.
- Changing to `seed_base = 1` changed both outputs.

Same-seed hashes:

- Partons: `c84c75e797ed10455998b8c0a1acbc33c33e832c5478301c40af3b58a3773b09`
- Hadrons: `d0e16461359ed653add08b13afaf57d65b5f33874977b3ab0762d39edd4b5145`

Different-seed hashes:

- Partons: `207f872072be70dc5b11a8def465d407a2f8ded64198f10b9f0666bbcb6bc9ae`
- Hadrons: `9d9f050c3c83385afe2522e61bd1dd9fd5f29ae214f6f2efac8c31e6f9d594a5`

## Current Scope and Limits

For the config-driven `main` path, the standard full averaged-hydro quenched workflow is now reproducible in practice under rerun from the beginning.

What this does not yet provide:

- direct replay of a single event by event number
- guaranteed documentation parity for the legacy trigger-oriented path
- protection against output-file appending if the same `output_base` is reused

## Notes

- Keep work on this branch unless explicitly told otherwise.
- Do not mix unrelated refactors into this branch.

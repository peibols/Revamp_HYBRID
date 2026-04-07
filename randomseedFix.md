# randomseedFix

This branch is reserved for the random-seed and reproducibility cleanup in the HYBRID revamp.

## Goal

Make the simulation reproducible in a controlled and documented way by removing mixed RNG usage and making all random streams explicit.

## Current Findings

- The shower-generation PYTHIA instance is explicitly seeded as `33 + njob`.
- The internal HYBRID RNG uses `mt19937_64` seeded as `1346 + njob`.
- A separate PYTHIA instance is used in `LundGenerator` for hadronization and is not explicitly reseeded in the current code path.
- The `ebe_hydro = 1` path still uses `rand()` in `HYBRID.cc`, which breaks the otherwise cleaner RNG scheme.
- The current IPSAT call path also passes the wrong quantity into the sampler: it draws a random integer first and then uses that integer as the sampler range.
- Reproducibility is reasonably good for the standard `ebe_hydro = 0` path when rerunning from the beginning with the same inputs and environment.
- Exact replay of a single event by event number is not currently supported because RNG state is not checkpointed per event.

## Branch Scope

- Introduce an explicit seed hierarchy.
- Remove remaining `rand()` usage.
- Explicitly seed all PYTHIA instances.
- Record effective seeds in output or startup logs.
- Add a minimal reproducibility smoke test.

## Proposed Seed Model

- `seed_base` from config
- `seed_shower = seed_base + 1`
- `seed_hybrid = seed_base + 2`
- `seed_lund = seed_base + 3`

This keeps the streams deterministic but distinct.

## Immediate Task List

- Add `seed_base` to the config interface.
- Thread explicit seeds into `TreeGenerator`.
- Thread explicit seeds into `LundGenerator`.
- Replace the `rand()` draw in the IPSAT path with the existing `numrand` engine.
- Fix the IPSAT sampling call to pass the full list size into the sampler.
- Print or write the effective seeds for each run.
- Add a two-run comparison test with identical seed input.

## Notes

- Keep work on this branch unless explicitly told otherwise.
- Do not mix unrelated refactors into this branch.

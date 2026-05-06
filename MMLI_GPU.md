# MMLI GPU Branch Notes

Branch: `main_moliere_lres_GPU_integration`

This branch starts GPU acceleration work from the validated MMLI physics branch while
keeping the production `main` executable CPU-identical unless a GPU path is explicitly
used.

## Implemented

- Added a standalone OpenCL backend for batched hydro interpolation:
  - `GpuHydroOpenCL.h`
  - `GpuHydroOpenCL.cc`
  - `gpu_hydro_benchmark.cc`
- Added read-only `HydroProfile` accessors needed to copy the hydro grid to GPU memory.
- Added `./compiler.sh gpu_hydro_benchmark`.
- Updated `compiler.sh` to default to the pinned PYTHIA 8.315 installation when present:
  `/data/yjlee/pythia/pythia8/pythia8315`

This first GPU step is intentionally limited to the hydro lookup primitive. It does not
alter event generation, energy-loss physics, Moliere scattering, LRES topology, wake, or
Lund hadronization.

## Local Benchmark

Host GPU:

- `NVIDIA GeForce GTX 1080 Ti`
- 11264 MiB

Hydro input:

- Averaged 0-5% hydro
- `/raid5/data/yjlee/hybrid_dev/hydro/hydro5020/hydroinfoPlaintxtHuichaoFormat_C0-5.dat`

Benchmark directory:

- `/raid5/data/yjlee/hybrid_dev/test/mmli_gpu_hydro_benchmark_20260506`

Commands:

```bash
./compiler.sh gpu_hydro_benchmark
mkdir -p /raid5/data/yjlee/hybrid_dev/test/mmli_gpu_hydro_benchmark_20260506
ln -sfn /raid5/data/yjlee/hybrid_dev/hydro/hydro5020/hydroinfoPlaintxtHuichaoFormat_C0-5.dat \
  /raid5/data/yjlee/hybrid_dev/test/mmli_gpu_hydro_benchmark_20260506/hydroinfoPlaintxtHuichaoFormat.dat
cd /raid5/data/yjlee/hybrid_dev/test/mmli_gpu_hydro_benchmark_20260506
/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_GPU_integration/gpu_hydro_benchmark --samples 10000 --repeats 10
/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_GPU_integration/gpu_hydro_benchmark --samples 100000 --repeats 3
/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_GPU_integration/gpu_hydro_benchmark --samples 1000000 --repeats 5
/raid5/data/yjlee/hybrid_dev/wt_main_moliere_lres_GPU_integration/gpu_hydro_benchmark --samples 5000000 --repeats 3
```

Results include host/device transfer for the query arrays and output arrays, with the
hydro grid resident on the GPU after initialization.

| Samples | CPU ms/repeat | GPU ms/repeat | Speedup | Max printed diff |
|---:|---:|---:|---:|---:|
| 10k | 0.908691 | 1.281352 | 0.709x | 0 |
| 100k | 17.783945 | 4.985805 | 3.567x | 0 |
| 1M | 175.701750 | 24.847081 | 7.071x | 0 |
| 5M | 933.325274 | 142.486390 | 6.550x | 0 |

Interpretation:

- Small batches are slower on GPU because transfer overhead dominates.
- Large batches show a conservative 6-7x speedup for hydro interpolation.
- This is not yet a full event-generation speedup. The next step is to batch `loss_rate`
  trajectory points or LRES segment propagation so full MMLI can actually use this kernel.

## Validation

- `./compiler.sh gpu_hydro_benchmark` builds.
- `./compiler.sh main` builds with pinned PYTHIA 8.315.
- The benchmark reports matching CPU/GPU values to the printed precision for temperature,
  `vx`, and `vy`.

## Next Integration Step

Wire the GPU backend into a guarded runtime path:

1. Add a config flag, for example `use_gpu_hydro = false`.
2. Keep the default path CPU-only and bitwise unchanged.
3. In `EnergyLoss::loss_rate`, collect hydro query points in batches and evaluate them via
   `GpuHydroOpenCL`.
4. Validate with `use_gpu_hydro=false` byte-identical to MMLI and `use_gpu_hydro=true`
   physics-equivalent in 1, 10, 100 event tests.

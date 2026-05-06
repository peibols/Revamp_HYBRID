#!/bin/bash

# Default to the pinned PYTHIA 8.315 installation used by validated MMLI campaigns.
# Override by setting PYTHIA_INCLUDE/PYTHIA_LIB before running.
PINNED_PYTHIA=${PINNED_PYTHIA:-/data/yjlee/pythia/pythia8/pythia8315}
if [[ -d "${PINNED_PYTHIA}/include" && -d "${PINNED_PYTHIA}/lib" ]]; then
  PYTHIA_INCLUDE=${PYTHIA_INCLUDE:-${PINNED_PYTHIA}/include}
  PYTHIA_LIB=${PYTHIA_LIB:-${PINNED_PYTHIA}/lib}
else
  PYTHIA_INCLUDE=${PYTHIA_INCLUDE:-/home/usc/ie/dpa/software/pythia-install/include}
  PYTHIA_LIB=${PYTHIA_LIB:-/home/usc/ie/dpa/software/pythia-install/lib}
fi

CXXFLAGS="-g -std=c++17 -Wall"

# Use xcrun on macOS to ensure proper SDK include paths (avoids libc++ header mismatch)
if [[ "$(uname)" == "Darwin" ]]; then
  CXX="$(xcrun --find clang++)"
  SYSROOT="$(xcrun --sdk macosx --show-sdk-path)"
  CXXFLAGS+=" -isysroot ${SYSROOT}"
else
  CXX="g++"
  CXXFLAGS+=" -mcmodel=medium"
fi

if [[ "$1" == "gpu_hydro_benchmark" ]]; then
  $CXX -O3 -std=c++17 -Wall \
    gpu_hydro_benchmark.cc HydroProfile.cc GpuHydroOpenCL.cc \
    -o gpu_hydro_benchmark \
    -ldl
  exit $?
fi

PYTHIA_INC_FLAGS=""
PYTHIA_LIB_FLAGS=""
GSL_LIB_FLAGS=""
if [[ -d "$PYTHIA_INCLUDE" ]]; then
  PYTHIA_INC_FLAGS="-I${PYTHIA_INCLUDE}"
fi
if [[ -d "$PYTHIA_LIB" ]]; then
  PYTHIA_LIB_FLAGS="-L${PYTHIA_LIB} -Wl,-rpath,${PYTHIA_LIB} -lpythia8 -ldl"
fi
if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists gsl; then
  GSL_LIB_FLAGS="$(pkg-config --libs gsl)"
else
  GSL_LIB_FLAGS="-lgsl -lgslcblas -lm"
fi

$CXX $CXXFLAGS \
	$1.cc TreeGenerator.cc Wake.cc Quench.cc Random.cc Parton.cc Hadron.cc \
	HydroProfile.cc WakeGenerator.cc LundGenerator.cc GlauberModel.cc EnergyLoss.cc HYBRID.cc Config.cc \
	MoliereTables.cc MoliereElastic.cc read_tables.cpp \
	-o $1 \
	$PYTHIA_INC_FLAGS $PYTHIA_LIB_FLAGS $GSL_LIB_FLAGS

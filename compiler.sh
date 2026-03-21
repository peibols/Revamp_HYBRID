#!/bin/bash

# Default PYTHIA paths (override by setting env vars PYTHIA_INCLUDE, PYTHIA_LIB before running).
PYTHIA_INCLUDE=${PYTHIA_INCLUDE:-/home/usc/ie/dpa/software/pythia-install/include}
PYTHIA_LIB=${PYTHIA_LIB:-/home/usc/ie/dpa/software/pythia-install/lib}

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

PYTHIA_INC_FLAGS=""
PYTHIA_LIB_FLAGS=""
if [[ -d "$PYTHIA_INCLUDE" ]]; then
  PYTHIA_INC_FLAGS="-I${PYTHIA_INCLUDE}"
fi
if [[ -d "$PYTHIA_LIB" ]]; then
  PYTHIA_LIB_FLAGS="-L${PYTHIA_LIB} -Wl,-rpath,${PYTHIA_LIB} -lpythia8 -ldl"
fi

$CXX $CXXFLAGS \
	$1.cc TreeGenerator.cc Wake.cc Quench.cc Random.cc Parton.cc Hadron.cc \
	HydroProfile.cc WakeGenerator.cc LundGenerator.cc GlauberModel.cc EnergyLoss.cc HYBRID.cc Config.cc \
	-o $1 \
	$PYTHIA_INC_FLAGS $PYTHIA_LIB_FLAGS

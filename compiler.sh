#!/bin/bash

PYTHIA_BIN=/home/usc/ie/dpa/software/pythia-install/bin
PYTHIA_INCLUDE=/home/usc/ie/dpa/software/pythia-install/include
PYTHIA_LIB=/home/usc/ie/dpa/software/pythia-install/lib
PYTHIA_SHARE=/home/usc/ie/dpa/software/pythia-install/share/Pythia8

g++ -g -std=c++0x -mcmodel=medium -Wall \
	$1.cc  Tree.cc HadWake.cc Wake.cc dEdx.cc SonsMomenta.cc Eloss.cc Quench.cc Glauber_PbPb.cc Hydro_PbPb.cc Random.cc Parton.cc Hadron.cc Lund.cc \
	Glauber_Chun.cc Hydro_Chun.cc \
	-o $1 \
	${PYTHIA_LIB}/libpythia8.a -I${PYTHIA_INCLUDE} -L${PYTHIA_LIB} -Wl,-rpath,${PYTHIA_LIB} -lpythia8 -ldl

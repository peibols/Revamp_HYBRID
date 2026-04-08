#!/bin/bash

g++ -g -Wall -std=c++11 -mcmodel=medium \
	$1.cc  read_tables.cpp Tree.cc HadWake.cc Wake.cc dEdx.cc SonsMomenta.cc Eloss.cc Quench.cc Glauber_PbPb.cc Hydro_PbPb.cc Random.cc Parton.cc Hadron.cc Lund.cc -o $1 \
	${PYTHIA8}/lib/libpythia8.a -I${PYTHIA8}/include -L${PYTHIA8}/lib -Wl,-rpath,${PYTHIA8}/lib -lpythia8 -ldl \
	-lgsl -lgslcblas -lm

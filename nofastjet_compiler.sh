#!/bin/bash

PYTHIA_BIN=/home/usc/ie/dpa/software/pythia-install/bin
PYTHIA_INCLUDE=/home/usc/ie/dpa/software/pythia-install/include
PYTHIA_LIB=/home/usc/ie/dpa/software/pythia-install/lib
PYTHIA_SHARE=/home/usc/ie/dpa/software/pythia-install/share/Pythia8

g++ -g -std=c++0x -mcmodel=medium -Wall \
        $1.cc -o $1 \
        ${PYTHIA_LIB}/libpythia8.a -I${PYTHIA_INCLUDE} -L${PYTHIA_LIB} -Wl,-rpath,${PYTHIA_LIB} -lpythia8 -ldl \

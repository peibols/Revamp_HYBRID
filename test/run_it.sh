#!/usr/bin/env bash

njob=1			#Seed for pythia events generation.
Nev=1000			#Number of events.
alpha=0.404		#Strength of strongly coupled energy loss.
kappa=15 		#Strength of gaussian broadening.
mode=0			#Selects energy loss rate: strongly coupled (0), radiative (1) and collisional (2). Just choose 0.
tmethod=1		#Selects value of pseudocritical temperature Tc = 145 MeV with tmethod=1. thmethod=0 selects Tc=170 MeV. Quenching happens up to that temperature. Just select tmethod=1
cent="0-5"		#Selects value of centrality. Since you only have 0-5% hdyro file, do not change this.
do_quench="true"	#Set it to false if you want to run in vacuum mode, i.e. no energy loss.

ebe_hydro=0

rm SOURCE.dat

ln -s /home/usc/ie/dpa/Collimator/cDGLAP_manyev/EbE_pbpb_qhat_samples/hydros_00-05/job-0/NcollList.dat NcollList.dat
ln -s /home/usc/ie/dpa/Collimator/cDGLAP_manyev/EbE_pbpb_qhat_samples/hydros_00-05/job-0/evolution_all_xyeta.dat evolution_all_xyeta.dat

time ./main $njob $Nev $cent $kappa $alpha $tmethod $do_quench $mode $ebe_hydro

#!/usr/bin/env bash

njob=1			#Seed for pythia events generation.
Nev=10			#Number of events.
alpha=0.404		#Strength of strongly coupled energy loss.
kappa=0 		#Strength of gaussian broadening.
mode=0			#Selects energy loss rate: strongly coupled (0), radiative (1) and collisional (2). Just choose 0.
tmethod=1		#Selects value of pseudocritical temperature Tc = 145 MeV with tmethod=1. thmethod=0 selects Tc=170 MeV. Quenching happens up to that temperature. Just select tmethod=1
cent="0-5"		#Selects value of centrality. Since you only have 0-5% hdyro file, do not change this.
do_quench="true"	#Set it to false if you want to run in vacuum mode, i.e. no energy loss.

trigger_pt=40.
trigger_eta=1.44
trigger_id=22

time ./main_trigger $njob $Nev $cent $kappa $alpha $tmethod $do_quench $mode $trigger_pt $trigger_eta $trigger_id

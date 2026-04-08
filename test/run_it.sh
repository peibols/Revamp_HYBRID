#!/bin/bash

N=4000
alpha=0.37
kappa=4.
mode=0
tmethod=1
hadro_type=1
cent="0-5"
njob=0

time ./main $njob $N $cent $kappa $alpha $tmethod true $mode $hadro_type

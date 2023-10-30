#!/bin/bash

eta_min=$1
eta_max=$2
eta_inc=$3

N=$4


for eta in `LANG=en_US seq $eta_min $eta_inc $eta_max`
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=16 thevenin_shear.jl 0.8 0.7 5e-4 1.0 $eta SINUSOIDAL 100.0 100 $N BAOAB 2.5 2.0 0.5 > nohup.out &
echo $eta
done

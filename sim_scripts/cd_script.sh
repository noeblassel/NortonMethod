#!/bin/bash

Nmin=$1
Nmax=$2
method=$3

for N in `seq $Nmin 10 $Nmax`
do

if [ $method == "norton" ]

then
    nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=4 norton_cd.jl 2.5 1.0 1e-3 1.0 0.7 10.0 100000 $N BAOAB > nohup.out &
else
    nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=4 thevenin_cd.jl 2.5 1.0 1e-3 1.0 10.0  10.0 100000 $N BAOAB > nohup.out &
fi
done
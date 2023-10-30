for eta in 1000.0 1500.0 2000.0 2500.0 3000.0 3500.0 4000.0 
do
nohup /libre/blasseln/julia-1.8.2/bin/julia --threads=8 thevenin_cd.jl 2.5 1.0 1e-3 1.0 $eta 10.0 1000 400 BAOAB > nohup.out &
done
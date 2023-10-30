using Plots

using Plots
using LinearAlgebra

methods=["linear","constant","sinusoidal"]
metamethods=["norton","thevenin"]
markershapes=[:xcross,:cross,:xcross]
colors=[:red,:green,:blue]
normalizations=[-4/π^2,2/π,1/2]

n_regr=10

ρ=0.7
γ=1.0
N=1000

L=cbrt(N/ρ)

colors=Dict("thevenin"=>:blue,"norton"=>:red)



norton_lines=readlines("norton_sinusoidal.dat")[2:end]
norton_dat=reduce(hcat,[parse.(Float64,split(l," ")) for l in norton_lines])

thevenin_lines =readlines("thevenin_sinusoidal.dat")[2:end]
thevenin_dat = reduce(hcat,[parse.(Float64,split(l," ")) for l in thevenin_lines])

η_min_lin = 0.0
η_max_lin = 1.0

η_min_nonlin = 1.0
η_max_nonlin = 1000.0

r,λ = norton_dat[[1,3],:]
η,R = thevenin_dat[[1,2],:]

r_lin = r[η_min_lin .< λ .< η_max_lin]
λ_lin = r[η_min_lin .< λ .< η_max_lin]

η_lin = η[ η_min_lin .< η .< η_max_lin]
R_lin = R[η_min_lin .< η .< η_max_lin]

r_nonlin = r[η_min_nonlin .< λ .< η_max_nonlin]
λ_nonlin = r[η_min_nonlin .< λ .< η_max_nonlin]

η_nonlin = η[ η_min_nonlin .< η .< η_max_nonlin]
R_nonlin = R[η_min_nonlin .< η .< η_max_nonlin]



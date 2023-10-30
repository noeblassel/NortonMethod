#!/bin/env julia

using Plots

thevenin_file,norton_file=ARGS

_,dat_norton...=map(split,readlines(norton_file))

Ns_norton= Float64[]
avs_norton=Float64[]

for (filename,avg,std,_,σ2,_)=dat_norton
    (r,Nx,Ny,Nz)=match(r"^norton_forcing_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",filename)
    r,avg,σ2,std=parse.(Float64,[r,avg,σ2,std])
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    push!(Ns_norton,Nx*Ny*Nz)
    push!(avs_norton,σ2) # delta method
end

_,dat_thevenin...=map(split,readlines(thevenin_file))

Ns_thevenin= Float64[]
avs_thevenin=Float64[]

for (filename,avg,std,_,σ2,_)=dat_thevenin
    (η,Nx,Ny,Nz)=match(r"^thevenin_response_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",filename)
    η,avg,σ2,std=parse.(Float64,[η,avg,σ2,std])
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    push!(Ns_thevenin,Nx*Ny*Nz)
    push!(avs_thevenin,σ2)
end


norton_perm = sortperm(Ns_norton)
thevenin_perm=sortperm(Ns_thevenin)

Ns_norton = Ns_norton[norton_perm]
avs_norton = avs_norton[norton_perm]

Ns_thevenin = Ns_thevenin[thevenin_perm]
avs_thevenin = avs_thevenin[thevenin_perm]

n_tot = length(Ns_thevenin)
n_regr = 8

lg_Ns_norton , lg_Ns_thevenin, lg_avs_norton,lg_avs_thevenin = map(x -> log.(x[n_tot-n_regr+1:n_tot]), [Ns_norton,Ns_thevenin,avs_norton,avs_thevenin])

lg_Ns_norton .-= first(lg_Ns_norton)
lg_Ns_thevenin .-= first(lg_Ns_thevenin)
lg_avs_norton .-= first(lg_avs_norton)
lg_avs_thevenin .-= first(lg_avs_thevenin)

av_rate_norton = inv(lg_Ns_norton'*lg_Ns_norton)*lg_Ns_norton'*lg_avs_norton
av_rate_thevenin = inv(lg_Ns_thevenin'*lg_Ns_thevenin)*lg_Ns_thevenin'*lg_avs_thevenin

reg_line_norton(n) =  avs_norton[n_tot-n_regr+1] * ((n/Ns_norton[n_tot-n_regr+1]) ^ av_rate_norton)
reg_line_thevenin(n) =  avs_thevenin[n_tot-n_regr+1] * ((n/Ns_norton[n_tot-n_regr+1]) ^ av_rate_thevenin)

scatter(Ns_norton,avs_norton,xlabel="N",ylabel="Asymptotic variance",markershape=:xcross,color=:red,label="Norton ⟨λ⟩ (exponent ≈ $(round(av_rate_norton,digits=3)))",xaxis=:log,yaxis=:log)
scatter!(Ns_thevenin,avs_thevenin,label="NEMD ⟨R⟩ (exponent ≈ $(round(av_rate_thevenin,digits=3)))",markershape=:cross,color=:blue)
plot!(reg_line_norton,label="",linestyle=:dot,color=:red)
plot!(reg_line_thevenin,label="",linestyle=:dot,color=:blue)
savefig("eq_straight_avs.pdf")



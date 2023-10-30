#!/bin/env julia

using Plots

thevenin_file,norton_file=ARGS

_,dat_norton...=map(split,readlines(norton_file))


rs_norton= Float64[]
λs_norton = Float64[]
avs_norton=Float64[]
vars_norton = Float64[]

fs = 12

for (filename,avg,std,_,σ2,_)=dat_norton
    (r,Nx,Ny,Nz)=match(r"^norton_forcing_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",filename)
    r,λ,σ2,std=parse.(Float64,[r,avg,σ2,std])
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    push!(rs_norton,r)
    push!(λs_norton,λ)
    push!(avs_norton,σ2*r^2/λ^4) # delta method
    push!(vars_norton,std^2)
end

_,dat_thevenin...=map(split,readlines(thevenin_file))

ηs_thevenin= Float64[]
Rs_thevenin=Float64[]
avs_thevenin=Float64[]
vars_thevenin=Float64[]

for (filename,avg,std,_,σ2,_)=dat_thevenin
    (η,Nx,Ny,Nz)=match(r"^thevenin_response_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",filename)
    η,R,σ2,std=parse.(Float64,[η,avg,σ2,std])
    push!(ηs_thevenin,η)
    push!(Rs_thevenin,R)
    push!(avs_thevenin,σ2/η^2)
    push!(vars_thevenin,std^2)
end


norton_perm = sortperm(rs_norton)
thevenin_perm=sortperm(ηs_thevenin)

rs_norton = rs_norton[norton_perm]
λs_norton = λs_norton[norton_perm]
avs_norton = avs_norton[norton_perm]
vars_norton = vars_norton[norton_perm]

ηs_thevenin = ηs_thevenin[thevenin_perm]
Rs_thevenin = Rs_thevenin[thevenin_perm]
avs_thevenin = avs_thevenin[thevenin_perm]
vars_thevenin = vars_thevenin[thevenin_perm]


scatter(λs_norton,rs_norton,xlabel="Forcing",ylabel="Normalized response",label="Norton",markershape=:xcross,color=:green,legend=:bottomright,xlims=(0,1.1maximum(ηs_thevenin)),ylims=(0,1.1maximum(rs_norton)),legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs)
scatter!(ηs_thevenin,Rs_thevenin,label="NEMD",markershape=:cross,color=:red)
n_tot = length(λs_norton)


savefig("response_plot.pdf")

n_regr = 4

lg_λs_norton , lg_ηs_thevenin, lg_avs_norton,lg_avs_thevenin = map(x -> log.(x[n_tot-n_regr+1:n_tot]), [λs_norton,ηs_thevenin,avs_norton,avs_thevenin])

lg_λs_norton .-= first(lg_λs_norton)
lg_ηs_thevenin .-= first(lg_ηs_thevenin)
lg_avs_norton .-= first(lg_avs_norton)
lg_avs_thevenin .-= first(lg_avs_thevenin)

av_rate_norton = inv(lg_λs_norton'*lg_λs_norton)*lg_λs_norton'*lg_avs_norton
av_rate_thevenin = inv(lg_ηs_thevenin'*lg_ηs_thevenin)*lg_ηs_thevenin'*lg_avs_thevenin

reg_line_norton(n) =  avs_norton[n_tot-n_regr+1] * ((n/λs_norton[n_tot-n_regr+1]) ^ av_rate_norton)
reg_line_thevenin(n) =  avs_thevenin[n_tot-n_regr+1] * ((n/ηs_thevenin[n_tot-n_regr+1]) ^ av_rate_thevenin)

scatter(λs_norton,avs_norton,xlabel="Forcing",ylabel="Asymptotic variance",markershape=:xcross,color=:red,label="Norton, (exponent ≈ $(round(av_rate_norton,digits=3)))",xaxis=:log,yaxis=:log,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs)
scatter!(ηs_thevenin,avs_thevenin,label="Thévenin, (exponent ≈ $(round(av_rate_thevenin,digits=3)))",markershape=:cross,color=:blue)
plot!(reg_line_norton,label="",linestyle=:dot,color=:red)
plot!(reg_line_thevenin,label="",linestyle=:dot,color=:blue)
savefig("neq_avs.pdf")

n_regr = 4  

lg_λs_norton , lg_ηs_thevenin, lg_vars_norton,lg_vars_thevenin = map(x -> log.(x[n_tot-n_regr+1:n_tot]), [λs_norton,ηs_thevenin,vars_norton,vars_thevenin])

lg_λs_norton .-= first(lg_λs_norton)
lg_ηs_thevenin .-= first(lg_ηs_thevenin)
lg_vars_norton .-= first(lg_vars_norton)
lg_vars_thevenin .-= first(lg_vars_thevenin)
av_rate_norton = inv(lg_λs_norton'*lg_λs_norton)*lg_λs_norton'*lg_vars_norton
av_rate_thevenin = inv(lg_ηs_thevenin'*lg_ηs_thevenin)*lg_ηs_thevenin'*lg_vars_thevenin

reg_line_norton(n) =  vars_norton[n_tot-n_regr+1] * ((n/λs_norton[n_tot-n_regr+1]) ^ av_rate_norton)
reg_line_thevenin(n) =  vars_thevenin[n_tot-n_regr+1] * ((n/ηs_thevenin[n_tot-n_regr+1]) ^ av_rate_thevenin)

scatter(λs_norton,vars_norton,xlabel="Forcing",ylabel="Variance",markershape=:xcross,color=:red,label="Norton, (exponent ≈ $(round(av_rate_norton,digits=3)))",legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs)#,xaxis=:log,yaxis=:log)
scatter!(ηs_thevenin,vars_thevenin,label="Thévenin, (exponent ≈ $(round(av_rate_thevenin,digits=3)))",markershape=:cross,color=:blue)
plot!(reg_line_norton,label="",linestyle=:dot,color=:red)
plot!(reg_line_thevenin,label="",linestyle=:dot,color=:blue)
savefig("neq_vars.pdf")


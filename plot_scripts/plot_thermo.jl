
#!/bin/env julia

using Plots

nemd_file,norton_file=ARGS

_,dat_norton...=map(split,readlines(norton_file))

fs = 12

ρ = 0.7
F1 = 0.5
γ = 1.0 

η,r = 0.3, 0.1

Ns_norton= Float64[]
λs_norton=Float64[]
avs_norton=Float64[]
vars_norton = Float64[]
Ls_norton = Float64[]
rs_norton = Float64[]

for (filename,avg,std,_,σ2,_)=dat_norton
    (r,Nx,Ny,Nz)=match(r"^norton_forcing_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",filename)
    r,avg,σ2,std=parse.(Float64,[r,avg,σ2,std])
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    push!(Ns_norton,Nx*Ny*Nz)
    push!(λs_norton,avg)
    push!(rs_norton,r)
    push!(avs_norton,σ2)
    push!(vars_norton,std^2)
    push!(Ls_norton,Ny*cbrt(inv(ρ)))
end

_,dat_nemd...=map(split,readlines(nemd_file))

Ns_nemd= Float64[]
Rs_nemd=Float64[]
avs_nemd=Float64[]
vars_nemd=Float64[]
Ls_nemd=Float64[]
ηs_nemd=Float64[]

for (filename,avg,std,_,σ2,_)=dat_nemd
    (η,Nx,Ny,Nz)=match(r"^thevenin_response_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",filename)
    η,avg,σ2,std=parse.(Float64,[η,avg,σ2,std])
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    push!(Ns_nemd,Nx*Ny*Nz)
    push!(Rs_nemd,avg)
    push!(ηs_nemd,η)
    push!(avs_nemd,σ2)
    push!(vars_nemd,std^2)
    push!(Ls_nemd,Ny*cbrt(inv(ρ)))
end


norton_perm = sortperm(Ns_norton)
nemd_perm=sortperm(Ns_nemd)

Ns_norton = Ns_norton[norton_perm]
λs_norton = λs_norton[norton_perm]
rs_norton = rs_norton[norton_perm]
avs_norton = avs_norton[norton_perm]
vars_norton = vars_norton[norton_perm]
Ls_norton = Ls_norton[norton_perm]


shear_norton = @. ρ*(F1 * λs_norton /rs_norton - γ)*(Ls_norton /2π)^2
asymptotic_vars_norton  = @. avs_norton * (ρ*F1/rs_norton * (Ls_norton/2π)^2)^2

Ns_nemd = Ns_nemd[nemd_perm]
Rs_nemd = Rs_nemd[nemd_perm]
ηs_nemd = ηs_nemd[nemd_perm]
avs_nemd = avs_nemd[nemd_perm]
vars_nemd = vars_nemd[nemd_perm]
Ls_nemd=Ls_nemd[nemd_perm]

shear_nemd = @. ρ*(F1 * ηs_nemd / Rs_nemd - γ)*(Ls_nemd /2π)^2
asymptotic_vars_nemd = @. avs_nemd * (ρ*F1 *ηs_nemd * (Ls_nemd/2π)^2)^2 / Rs_nemd^4


#= println("N norton: ",Ns_norton)
println("ρ norton: ",us_norton)
println("σ2 norton: ",avs_norton)
println()
println("N nemd: ",Ns_nemd)
println("ρ nemd: ",us_nemd)
println("σ2 nemd: ",avs_nemd) =#

us_norton = rs_norton ./ λs_norton
us_nemd = Rs_nemd ./ ηs_nemd

scatter(Ns_norton,F1*inv.(us_norton).-γ,xlabel="N",ylabel="F₁/U₁-γₓ",label="Norton",markershape=:xcross,color=:green,legend=:topright,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs)#,yaxis=:log,xaxis=:log)
scatter!(Ns_nemd,F1*inv.(us_nemd).-γ,label="NEMD",markershape=:cross,color=:red)

### least squares extrapolation to large N limit ###
n_tot = length(Ns_norton)
n_regr = 8

lg_Ns_norton =log.(Ns_norton[n_tot-n_regr+1:n_tot])
X = ones(n_regr,2)
X[:,2] .= lg_Ns_norton
Y = log.(F1*inv.(us_norton[n_tot-n_regr+1:n_tot]).-γ)

Theta = inv(X'*X)*(X'*Y)
a,b = Theta

println(Theta, X\Y)

println(a,b)
println("norton η ≈ ",a * ρ^(2/3)/π^2)

#hline!([a], label="Extrapolated limit ($(round(a,digits=3)))",color=:black)
plot!(t->exp(a)*t^b,label="",color=:green,linestyle=:dash)#,ylims=(0.9*minimum(us_norton),1.1*maximum(us_norton)))

lg_Ns_nemd = log.(Ns_nemd[n_tot-n_regr+1:n_tot])
X = ones(n_regr,2)
X[:,2] .= lg_Ns_nemd
Y = log.(F1*inv.(us_nemd[n_tot-n_regr+1:n_tot]).-γ)


c,d = inv(X'*X)*(X'*Y)
println(c,d)
println("nemd ",d * ρ^(2/3)/π^2)
plot!(t->exp(c)*t^d,label="",color=:red,linestyle=:dot)

savefig("neq_thermo_limit.pdf")

scatter(Ns_norton,F1*inv.(us_norton).-γ,xlabel="N",ylabel="F₁/U₁-γₓ",label="Norton (slope = $(round(b,digits=3)))",markershape=:xcross,color=:green,legend=:topright,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs,yaxis=:log,xaxis=:log)
scatter!(Ns_nemd,F1*inv.(us_nemd).-γ,label="NEMD (slope = $(round(d,digits=3)))",markershape=:cross,color=:red)
plot!(t->exp(a)*t^b,label="",color=:green,linestyle=:dash)
plot!(t->exp(c)*t^d,label="",color=:red,linestyle=:dot)

savefig("neq_thermo_limit_log.pdf")

n_regr = 4

lg_Ns_norton , lg_Ns_nemd, lg_avs_norton,lg_avs_nemd = map(x -> log.(x[n_tot-n_regr+1:n_tot]), [Ns_norton,Ns_nemd,asymptotic_vars_norton,asymptotic_vars_nemd])

lg_Ns_norton .-= first(lg_Ns_norton)
lg_Ns_nemd .-= first(lg_Ns_nemd)
lg_avs_norton .-= first(lg_avs_norton)
lg_avs_nemd .-= first(lg_avs_nemd)

av_rate_norton = inv(lg_Ns_norton'*lg_Ns_norton)*lg_Ns_norton'*lg_avs_norton
av_rate_nemd = inv(lg_Ns_nemd'*lg_Ns_nemd)*lg_Ns_nemd'*lg_avs_nemd

reg_line_norton(n) =  avs_norton[n_tot-n_regr+1] * ((n/Ns_norton[n_tot-n_regr+1]) ^ av_rate_norton)
reg_line_nemd(n) =  avs_nemd[n_tot-n_regr+1] * ((n/Ns_norton[n_tot-n_regr+1]) ^ av_rate_nemd)

scatter(Ns_norton,asymptotic_vars_norton,xlabel="N",ylabel="Asymptotic variance",markershape=:xcross,color=:red,label="Norton, r = $r, (asymptotic exponent ≈ $(round(av_rate_norton,digits=3)))",xaxis=:log,yaxis=:log,legend=:right,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs)
scatter!(Ns_nemd,asymptotic_vars_nemd,label="NEMD, η = $η, (asymptotic exponent ≈ $(round(av_rate_nemd,digits=3)))",markershape=:cross,color=:blue)
#plot!(reg_line_norton,label="",linestyle=:dot,color=:red)
#plot!(reg_line_nemd,label="",linestyle=:dot,color=:blue)
savefig("neq_avs.pdf")

n_regr = 7

lg_Ns_norton , lg_Ns_nemd, lg_vars_norton,lg_vars_nemd = map(x -> log.(x[n_tot-n_regr+1:n_tot]), [Ns_norton,Ns_nemd,vars_norton,vars_nemd])

lg_Ns_norton .-= first(lg_Ns_norton)
lg_Ns_nemd .-= first(lg_Ns_nemd)
lg_vars_norton .-= first(lg_vars_norton)
lg_vars_nemd .-= first(lg_vars_nemd)

var_rate_norton = inv(lg_Ns_norton'*lg_Ns_norton)*lg_Ns_norton'*lg_vars_norton
var_rate_nemd = inv(lg_Ns_nemd'*lg_Ns_nemd)*lg_Ns_nemd'*lg_vars_nemd

reg_line_norton(n) =  vars_norton[n_tot-n_regr+1] * ((n/Ns_norton[n_tot-n_regr+1]) ^ var_rate_norton)
reg_line_nemd(n) =  vars_nemd[n_tot-n_regr+1] * ((n/Ns_norton[n_tot-n_regr+1]) ^ var_rate_nemd)

scatter(Ns_norton,vars_norton,xlabel="N",ylabel="Variance",markershape=:xcross,color=:red,label="Norton, r = $r, (exponent ≈ $(round(var_rate_norton,digits=3)))",xaxis=:log,yaxis=:log,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs)
scatter!(Ns_nemd,vars_nemd,label="NEMD, η = $η, (exponent ≈ $(round(var_rate_nemd,digits=3)))",markershape=:cross,color=:blue)
plot!(reg_line_norton,label="",linestyle=:dot,color=:red)
plot!(reg_line_nemd,label="",linestyle=:dot,color=:blue)
savefig("neq_vars.pdf")

scatter(Ns_norton,shear_norton)
scatter!(Ns_nemd,shear_nemd)
savefig("shear_thermo.pdf")

println("σ²(η) norton :",asymptotic_vars_norton)
println("σ²(η) nemd: ",asymptotic_vars_nemd)

println("σ²(λ) :", avs_norton)
println("σ²(R) :",avs_nemd)
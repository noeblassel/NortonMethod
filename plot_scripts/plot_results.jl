using Plots
using LinearAlgebra

methods=["linear","constant","sinusoidal"]
markershapes=[:xcross,:cross,:rect]
colors=[:red,:green,:blue]
normalizations=[-1,1,1]

n_regr=5

fs = 12
#legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs

lin_plots=Dict(method=>plot(legend=:bottomright,xlabel="Forcing",ylabel="Normalized response U₁/F₁",labelfontsize=fs,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs) for method in methods)
full_plots =Dict(method=>plot(legend=:bottomright,xlabel="Forcing",ylabel="Normalized response U₁/F₁",labelfontsize=fs,legendfontsize=fs,tickfontsize=fs,titlefontsize=fs) for method in methods)


ρ=0.7
γ=1.0

N=1000
L=cbrt(N/ρ)

fourier_coeffs = Dict("linear" => 4/π^2, "constant" => 2/π , "sinusoidal" => 1/2)

n_lin = 10

for (i,method) in enumerate(methods)
    println(method)
    lines=readlines("norton_$(method).dat")[2:end]
    dat=reduce(hcat,[parse.(Float64,split(l," ")) for l in lines])
    forcings=dat[3,:]
    responses=dat[1,:]

    if method=="linear"
        forcings .= abs.(forcings)
        responses .= abs.(responses)
    end

    show(forcings)
    show(responses)

    forcing_perm = sortperm(forcings)
    forcings = forcings[forcing_perm]
    responses = responses[forcing_perm]

    ρ_norton=inv(dot(forcings[1:n_regr],forcings[1:n_regr]))*dot(forcings[1:n_regr],responses[1:n_regr])
    println("\tTransport coefficient: ",ρ_norton)
    σ=ρ*(inv(ρ_norton)*fourier_coeffs[method]-γ)*(L/2π)^2
    println("\tShear viscosity :",σ)
    sigmas= sqrt.(dat[5,:])

    shears = @. ρ*((forcings * fourier_coeffs[method] / responses) - γ) * (L/2π)^2

    scatter!(full_plots[method],forcings[n_lin+1:end],responses[n_lin+1:end]/fourier_coeffs[method],label="Norton",color=:green,markershape=:xcross,markersize=8)
    scatter!(lin_plots[method],forcings,responses/fourier_coeffs[method],label="Norton ($(round(ρ_norton/fourier_coeffs[method],digits=3)))",color=:green,markershape=:xcross,markersize=8)

    xlims!(lin_plots[method],0,1.1*forcings[n_lin])
    ylims!(lin_plots[method],0,1.1*responses[n_lin]/fourier_coeffs[method])

    xlims!(full_plots[method],0,1.1*last(forcings))
    ylims!(full_plots[method],0,1.1*last(responses)/fourier_coeffs[method])

    hline(lin_plots[method],[σ],label="",linestyle=:dash,color=:green)

    lines=readlines("thevenin_$(method).dat")[2:end]
    dat=reduce(hcat,[parse.(Float64,split(l," ")) for l in lines])
    forcings=dat[1,:]

    if method=="linear"
        responses = abs.(dat[3,:])
    else
        responses = dat[2,:]
    end

    forcing_perm = sortperm(forcings)
    forcings = forcings[forcing_perm]
    responses = responses[forcing_perm]

    ρ_nemd=inv(dot(forcings[1:n_regr],forcings[1:n_regr]))*dot(forcings[1:n_regr],responses[1:n_regr])
    σ=ρ*(inv(ρ_nemd)*fourier_coeffs[method]-γ)*(L/2π)^2
    shears = @. ρ*((forcings * fourier_coeffs[method] / responses) - γ) * (L/2π)^2

    scatter!(full_plots[method],forcings[n_lin+1:end],responses[n_lin+1:end]/fourier_coeffs[method],label="NEMD",color=:red,markershape=:cross,markersize=8)
    scatter!(lin_plots[method],forcings,responses/fourier_coeffs[method],label="NEMD ($(round(ρ_nemd/fourier_coeffs[method],digits=3)))",color=:red,markershape=:cross,markersize=8)


    

    #hline!(lin_plots[method],[σ],label="",linestyle=:dot,color=:red)
    plot!(lin_plots[method],label="fitted linear response",color=:black,linestyle=:dot, t->ρ_norton*t/fourier_coeffs[method])
    plot!(full_plots[method],label="fitted linear response",color=:black,linestyle=:dot, t->ρ_norton*t/fourier_coeffs[method])
    savefig(lin_plots[method],"shears_lin_$(method).pdf")
    savefig(full_plots[method],"shears_nonlin_$(method).pdf")
end

    

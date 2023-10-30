using Plots
using LinearAlgebra

metamethods=["norton","thevenin"]
methods=["linear","constant","sinusoidal"]
markershapes=[:xcross,:cross,:xcross]
colors=[:red,:green,:blue]
normalizations=[-4/π^2,2/π,1/2]

n_regr=10



ρ=0.7
γ=1.0

N=1000
L=cbrt(N/ρ)

fs = 12
full_plot=plot(legend=:topright,xlabel="Forcing magnitude",ylabel=" Asymptotic variance for U₁",xaxis=:log,yaxis=:log, markerfontsize=fs,labelfontsize=fs,tickfontsize=fs)

forcing_ixs=Dict(("norton","sinusoidal")=>3,("norton","constant")=>3,("norton","linear")=>3,("thevenin","sinusoidal")=>1,("thevenin","constant")=>1,("thevenin","linear")=>1)
response_ixs=Dict(("norton","sinusoidal")=>1,("norton","constant")=>1,("norton","linear")=>1,("thevenin","sinusoidal")=>2,("thevenin","constant")=>2,("thevenin","linear")=>3)
avs_ixs=Dict(("norton","sinusoidal")=>5,("norton","constant")=>5,("norton","linear")=>5,("thevenin","sinusoidal")=>5,("thevenin","constant")=>5,("thevenin","linear")=>6)

for metamethod in metamethods
    println("$metamethod")
    shape = (metamethod =="norton") ? :xcross : :cross
    for (i,method) in enumerate(methods)
        println(method)
        lines=readlines("$(metamethod)_$(method).dat")[2:end]
        dat=reduce(hcat,[parse.(Float64,split(l," ")) for l in lines])

        forcings=dat[forcing_ixs[(metamethod,method)],:]
        avs=dat[avs_ixs[(metamethod,method)],:]
        Rs=dat[response_ixs[(metamethod,method)],:]

        if method=="linear"
            forcings .= abs.(forcings)
            Rs .= abs.(Rs)
        end

        if metamethod == "norton"
           @. avs = ((Rs^ 2) * avs / forcings^4) #delta method
        else
            @. avs = avs / forcings^2
        end

        name = (metamethod == "norton") ? "Norton" : "NEMD"

        scatter!(full_plot,forcings,avs,label="$(name) ($(method) force profile)",markershape=shape)

    end
end

plot!(full_plot,t-> 0.1inv(t^2),linestyle=:dot, color=:red,label="slope=-2")
savefig(full_plot,"norton_full_avs.pdf")
    

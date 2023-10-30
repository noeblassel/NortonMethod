using Plots,Statistics

nemd_file,norton_file = ARGS

λs = Float64[]
σλs = Float64[]
avλs = Float64[]
Nλs = Float64[]
rλs = Float64[]

Rs = Float64[]
σRs = Float64[]
avRs = Float64[]
NRs=Float64[]
ηRs = Float64[]

lines_nemd = readlines(nemd_file)[2:end]
lines_norton = readlines(norton_file)[2:end]

for l in lines_nemd
    f,m,σ,_,av,_ = split(l)
    (m,σ,av) = parse.(Float64,[m,σ,av])
    push!.([Rs,σRs,avRs],[m,σ,av])

    N,ρ,T,η =  parse.(Float64, match(r"nemd_response_colordrift_(.+)_(.+)_(.+)_(.+).out",f))
    push!.([NRs,ηRs],[N,η])

end

for l in lines_norton
    f,m,σ,_,av,_ = split(l)
    (m,σ,av) = parse.(Float64,[m,σ,av])
    push!.([λs,σλs,avλs],[m,σ,av])

    N,ρ,T,r =  parse.(Float64, match(r"norton_forcing_colordrift_(.+)_(.+)_(.+)_(.+).out",f))
    push!.([Nλs,rλs],[N,r])
end

nemd_perm = sortperm(NRs)
norton_perm = sortperm(Nλs)

λs=λs[norton_perm]
σλs=σλs[norton_perm]
avλs=avλs[norton_perm]
Nλs=Nλs[norton_perm]
rλs = rλs[norton_perm]

Rs=Rs[nemd_perm]
σRs=σRs[nemd_perm]
avRs=avRs[nemd_perm]
NRs=NRs[nemd_perm]
ηRs = ηRs[nemd_perm]

resp_plot = plot(xlabel="N",ylabel="estimated mobility")
scatter!(resp_plot,NRs,Rs ./ ηRs, markershape=:xcross, color=:red, label="nemd")
scatter!(resp_plot,Nλs,rλs ./ λs, markershape=:cross, color=:green, label="norton")

Xnemd = ones(length(NRs),2)
Xnemd[:,2] .= inv.(NRs)

Ynemd = Rs ./ ηRs

a_nemd, b_nemd = inv(Xnemd'*Xnemd)*(Xnemd'*Ynemd)

plot!(resp_plot,n-> a_nemd + b_nemd/n,label="",color=:black,linestyle=:dot)

Xnorton = ones(length(Nλs),2)
Xnorton[:,2] .= inv.(Nλs)

Ynorton = rλs ./ λs

a_norton, b_norton = inv(Xnorton'*Xnorton)*(Xnorton'*Ynorton)
plot!(resp_plot,n-> a_norton + b_norton/n,label="",color=:black,linestyle=:dot)


println("nemd : $(a_nemd) +$(b_nemd)/N ; norton : $(a_norton)+$(b_norton)/N")

var_plot = plot(xlabel="N",ylabel="variance of estimator",xaxis=:log,yaxis=:log)
scatter!(var_plot,NRs,σRs, markershape=:xcross,color=:red,label="nemd")
scatter!(var_plot,Nλs,σλs, markershape=:cross,color=:green,label="norton")




savefig(resp_plot,"thermo_cd.pdf")
savefig(var_plot,"vars_cd.pdf")


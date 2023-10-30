using Plots,LinearAlgebra

include("norton_vars.out")
Ns_norton  = Ns
vars_norton = vars

include("thevenin_vars.out")
Ns_thevenin = Ns
vars_thevenin=vars

pl=plot(xlabel="N",ylabel="variance",xaxis=:log,yaxis=:log)

vars_norton = vars_norton[sortperm(Ns_norton)]
sort!(Ns_norton)

lg_Ns_norton = log.(Ns_norton) .- log(Ns_norton[1])
lg_vars_norton = log.(vars_norton) .- log(vars_norton[1])

lg_Ns_thevenin = log.(Ns_thevenin) .- log(Ns_thevenin[1])
lg_vars_thevenin = log.(vars_thevenin) .- log(vars_thevenin[1])

regression_coeff_norton = inv(dot(lg_Ns_norton,lg_Ns_norton)) * dot(lg_Ns_norton,lg_vars_norton)
regression_coeff_thevenin = inv(dot(lg_Ns_thevenin,lg_Ns_thevenin))*dot(lg_Ns_thevenin,lg_vars_thevenin)

scatter!(pl,Ns_norton,vars_norton,label="λ (slope=$(round(regression_coeff_norton,digits=2))",color=:blue,markershape=:xcross)
scatter!(pl,Ns_thevenin,vars_thevenin,label="R (slope=$(round(regression_coeff_thevenin,digits=2))",color=:red,markershape=:xcross)

plot!(pl,n->vars_norton[1]*(n/Ns_norton[1])^regression_coeff_norton,label="",color=:blue,linestyle=:dot)
plot!(pl,n->vars_thevenin[1]*(n/Ns_thevenin[1])^regression_coeff_thevenin,label="",color=:red,linestyle=:dot)

savefig(pl,"vars_zero_ensemble.pdf")

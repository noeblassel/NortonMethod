#!/bin/env julia
using Plots


println("Plots two histograms. Usage: ./plot_two_hists.jl HISTOGRAM_FILE1 HISTOGRAM_FILE2 LABEL1 LABEL2 XLABEL YLABEL")
println("Histogram files should be structured as follows:")
println("Line 1: min_value, Max_value, N_bins")
println("Line 2: N_bins comma-separated values")

hist1,hist2,lbl1,lbl2,xlabel,ylabel=ARGS

pl=plot(xlabel=xlabel,ylabel=ylabel)

for (hist,lbl) in zip([hist1,hist2],[lbl1,lbl2])
    a,b=readlines(hist)
    m,M,N=split(a)

    m,M=parse.(Float64,(m,M))
    N=parse(Int64,N)

    R=range(m,M,N)
    dx=R.step

    bins=parse.(Float64,split(b,","))
    S=sum(bins)
    plot!(pl,R,bins/Float64(S*dx),label=lbl)
end

savefig("joint_$(hist1)_$(hist2).pdf")

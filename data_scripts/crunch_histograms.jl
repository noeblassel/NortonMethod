#!/bin/env julia

using Plots

output_file,histogram_files... = ARGS

println("usage: ./crunch_histograms.jl output_file histogram_files")

vars = Float64[]
Ns = Int64[]
file_regex=r"histogram_1000_thevenin_response_SINUSOIDAL_0\.0_(\d+)_(\d+)_(\d+).out"
for file in histogram_files
    println("crunching $file")

    a,b,c=parse.(Int64,match(file_regex,file))
    push!(Ns,a*b*c)

    lines=readlines(file)
    l,h,n_bins=split(lines[1])
    hist=parse.(Float64,split(lines[2],","))
    l=parse(Float64,l)
    h=parse(Float64,h)
    n_bins=parse(Int64,n_bins)
    dx=(h-l)/n_bins
    n_samps=sum(hist)
    sum_hist = sum_sq_hist = 0
    for (i,v)=enumerate(range(l+dx/2,h-dx/2,n_bins))
        sum_hist += v*hist[i]
        sum_sq_hist += v^2*hist[i]
    end
    σ2 = sum_sq_hist/n_samps - (sum_hist/n_samps)^2
    push!(vars,σ2)
end

f= open(output_file,"w")
println(f,"Ns=",Ns)
println(f,"vars=",vars)
close(f)
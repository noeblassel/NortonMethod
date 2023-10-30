
method = ARGS[1]
r= ARGS[2]
using Plots

Ns=[3,4,5,6,7,8,9,10]
n_sigma = 3

fs=12

file_prefix="histogram_1000_$(method=="norton" ? "norton_forcing" : "thevenin_response")_SINUSOIDAL_$(r)_"
file_suffix=".out"

xlims = (method=="norton") ?  (-0.2,0.2) : (-0.1,0.1)

means=zero(Ns)
mean_sqs=zero(Ns)

xlabel = (method == "norton") ? "λ" : "R"
hist_plot=plot(xlabel=xlabel,ylabel="probability density",legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs,xlims=xlims)
std_hist_plot=plot(xlabel=xlabel,ylabel="standardized probability density",legendfontsize=fs,tickfontsize=fs,titlefontsize=fs,labelfontsize=fs,xlims=xlims)

std_range=range(-5,5,1000)


min_low = Inf
max_high = -Inf

for N in Ns
    println(N)
    file=open(string(file_prefix,"$(2N)_$(4N)_$(N)",file_suffix),"r")
    lines=readlines(file)
    l,h,n_bins,avg,σ=split(lines[1])
    hist=parse.(Float64,split(lines[2],","))
    l=parse(Float64,l)
    h=parse(Float64,h)
    global min_low=min(l,min_low)
    global max_high=max(h,max_high)
    n_bins=parse(Int64,n_bins)
    dx=(h-l)/n_bins

    dx_norm = 10/1000
    hist_std = hist / (dx_norm * sum(hist))

    hist /= dx*sum(hist)
    r=range(l,h,n_bins)
    plot!(hist_plot,r,hist,label="N=$(8N^3)")
    plot!(std_hist_plot,std_range,hist_std,label="N=$(8N^3)")
end


max_10σ=max_high-min_low
center= (max_high+min_low)/2
max_nσ=max_10σ*n_sigma/10

#xlims!(hist_plot, (center-max_nσ,center+max_nσ))

savefig(hist_plot,"histograms_neq_$(method)_$(r).pdf")
savefig(std_hist_plot,"std_histogram_neq_$(method)_$(r).pdf")
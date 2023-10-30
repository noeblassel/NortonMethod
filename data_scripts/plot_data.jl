
using Plots

norton_dir="./norton"
thevenin_dir="./thevenin"
hist_dirname="histograms"

hists_norton=readdir(joinpath(norton_dir,hist_dirname))
hists_thevenin=readdir(joinpath(thevenin_dir,hist_dirname))

output_file=open("moments.txt","w")

Ns=[512,1000,1728,2744]
color_palette=palette(:default)

function moments(m,M,bins)
    N_bins=length(bins)
    xrange=range(m,M,N_bins)
    S=sum(bins)
    avg=sum(x1*b1+x2*b2 for (x1,x2,b1,b2)=zip(xrange[1:end-1],xrange[2:end],bins[1:end-1],bins[2:end]))/2S
    var=sum((x1-avg)^2*b1+(x2-avg)^2*b2 for (x1,x2,b1,b2)=zip(xrange[1:end-1],xrange[2:end],bins[1:end-1],bins[2:end]))/2S
    return avg,var
end

println(output_file,"NORTON RESULTS: N r λ var_λ")

norton_plot=plot(xlabel="λ",ylabel="likelihood")
av_plot=plot(xlabel="N",ylabel="Asymptotic variance",xaxis=:log,yaxis=:log)

for hist in hists_norton
    r,Nx,Ny,Nz=match(r"^histogram_1000_norton_forcing_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",hist)
    r=parse(Float64,r)
    Nx,Ny,Nz=parse.(Int,[Nx,Ny,Nz])
    println(hist)
    f=open(joinpath(norton_dir,hist_dirname,hist))
    m,M,Nbins=split(readline(f))
    m,M=parse.(Float64,[m,M])
    Nbins=parse(Int,Nbins)
    bins=parse.(Float64,split(readline(f),','))
    λ,var_λ=moments(m,M,bins)
    N=Nx*Ny*Nz
    println(output_file,"$N $r $λ $var_λ")
    
    if N in Ns
        i=findfirst(N .== Ns)
        println(i)
        xrange=range(m,M,Nbins)
        dx=xrange.step.hi
        S=sum(bins)
        plot!(norton_plot,xrange,bins/(S*dx),label="$N",color=color_palette[i])
        vline!([λ],linestyle=:dot,label="",color=color_palette[i])
    end
end


thevenin_plot=plot(xlabel="R",ylabel="likelihood")

savefig(norton_plot,"norton_plot.pdf")

println(output_file,"THEVENIN RESULTS: N η R var_R")

for hist in hists_thevenin
    η,Nx,Ny,Nz=match(r"^histogram_1000_thevenin_response_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",hist)
    η=parse(Float64,η)
    Nx,Ny,Nz=parse.(Int,[Nx,Ny,Nz])
    println(hist)
    f=open(joinpath(thevenin_dir,hist_dirname,hist))
    m,M,Nbins=split(readline(f))
    m,M=parse.(Float64,[m,M])
    Nbins=parse(Int,Nbins)
    bins=parse.(Float64,split(readline(f),','))
    R,var_R=moments(m,M,bins)
    N=Nx*Ny*Nz
    println(output_file,"$N $η $R $var_R")
    if N in Ns
        i=findfirst(N .== Ns)
        println(i)
        xrange=range(m,M,Nbins)
        dx=xrange.step.hi
        S=sum(bins)
        plot!(thevenin_plot,xrange,bins/(S*dx),label="$N",color=color_palette[i])
        vline!([R],label="",linestyle=:dot,color=color_palette[i])
    end
end

savefig(thevenin_plot,"thevenin_plot.pdf")
close(output_file)

_,lines_norton... = readlines(joinpath(norton_dir,"forcing_results.txt"))
_,lines_thevenin... = readlines(joinpath(thevenin_dir,"response_results.txt"))

avs_norton=Float64[]
avs_thevenin=Float64[]
Ns_norton=Float64[]
Ns_thevenin=Float64[]

for l in lines_norton
    pathname,m,_,σ2,_=split(l)
    r,Nx,Ny,Nz=match(r"^.+_(.+)_(\d+)_(\d+)_(\d+).out$",pathname)
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    N=Nx*Ny*Nz
    r,m,σ2=parse.(Float64,[r,m,σ2])
    av_ρ = σ2*r^2/m^4
    push!(avs_norton,av_ρ)
    push!(Ns_norton,N)
end

for l in lines_thevenin
    pathname,m,_,σ2,_=split(l)
    η,Nx,Ny,Nz=match(r"^.+_(.+)_(\d+)_(\d+)_(\d+).out$",pathname)
    Nx,Ny,Nz=parse.(Int64,[Nx,Ny,Nz])
    N=Nx*Ny*Nz
    η,m,σ2=parse.(Float64,[η,m,σ2])
    av_ρ = σ2/η^2
    push!(avs_thevenin,av_ρ)
    push!(Ns_thevenin,N)
end


logNt=log.(Ns_thevenin)
logavt=log.(avs_thevenin)

logNt .-= first(logNt)
logavt .-= first(logavt)

α_t=inv(logNt'*logNt) * (logNt' * logavt)
b_t = first(avs_thevenin)/first(Ns_thevenin)^α_t

scatter!(av_plot,Ns_norton,avs_norton,markershape=:xcross,label="Norton", color=:red)
scatter!(av_plot,Ns_thevenin,avs_thevenin,markershape=:xcross,label="Thévenin (slope: $(round(α_t,digits=2)))",color=:blue)
plot!(av_plot,x->b_t*x^α_t,linestyle=:dot,color=:blue,label="")
savefig(av_plot,"asymptotic_variances.pdf")
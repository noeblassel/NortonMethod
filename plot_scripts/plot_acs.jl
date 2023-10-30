using Plots

file = ARGS[1]
regex = Regex(ARGS[2])
ncorr = parse(Int64,ARGS[3])

pl=plot(xlabel="time",ylabel="autocorrelation")
dt=1e-3

trange = dt*(0:ncorr-1)

for l in readlines(file)[2:end]
    fname,m,msq,series=split(l)
    Nx,Ny,Nz=parse.(Int64,match(regex,fname))
    N=Nx*Ny*Nz
    m,msq=parse.(Float64,[m,msq])
    series=parse.(Float64,split(series,","))
    plot!(pl,trange,series[1:ncorr],label="N=$N")
end

savefig(pl,file*".pdf")
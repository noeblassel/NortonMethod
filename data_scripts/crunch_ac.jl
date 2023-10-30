using Statistics, Base.Threads

n_corr, output_file, base_dir ,files... = ARGS

n_corr = parse(Int64,n_corr)

function crunch_autocorrelation(filename,n_corr)
    series = reinterpret(Float64,read(filename))
    m=mean(series)
    series .-= m
    msq = mean(series .^ 2)

    C_threaded=[zeros(n_corr) for i=1:nthreads()]
    n_done=[0 for i=1:nthreads()]

    @threads for i=n_corr+1:length(series)
        C_threaded[threadid()] +=  series[i] * series[i:-1:i-n_corr+1]
        n_done[threadid()]+=1

        if sum(n_done) % 10000 == 0
            println("$(sum(n_done)) / $(length(series))")
        end
    end

    C=sum(C_threaded)
    C/= (length(series) - n_corr -1)

    return C, m, msq
end


of=open(output_file,"w")

println(of,"FILENAME MEAN MEAN_SQUARE AUTOCORRELATION_SERIES")

for (i,f) in enumerate(files)
    println("$i / $(length(files))")
    C,m,msq = crunch_autocorrelation(joinpath(base_dir,f), n_corr)
    println(of,"$f $m $msq $(join(C,','))")
end

close(of)
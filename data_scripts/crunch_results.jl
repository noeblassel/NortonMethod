using Statistics

println("USAGE: N_BINS_HISTOGRAM N_σ OUTPUT_FILENAME INPUT_FILES[]")
n_bins,n_σ, output_file,input_files... = ARGS
n_bins=parse(Int64,n_bins)
n_σ=parse(Float64,n_σ)

if !isdir("histograms")
    mkdir("histograms")
end

function asymptotic_variance(v, min_n_blocks = 10)
    data=copy(v)
    avg=mean(v)
    data .-= avg

    if length(data)%2 ==1
        pop!(data)
    end

    N=length(data)
    K = 1 # block size
    max_var=0.0

    while length(data) > min_n_blocks
        var_K = var(data) * K

        if var_K > max_var
            max_var = var_K
        else
            break 
        end

        K *=2
        data = (data[1:2:end]+data[2:2:end])/2

        if length(data)%2 ==1 # always keep an even number of blocks
            pop!(data)
        end

    end
    return max_var
end

header = "FILENAME MEAN STD N_SAMPLES ASYMPTOTIC_VARIANCE HISTOGRAM_PATH"
f_output=open(output_file,"w")

println(f_output,header)

for f in input_files
    println("Processing $f:")

    try
        series=reinterpret(Float64,read(f))
        av=asymptotic_variance(series)
        avg=mean(series)
        σ = std(series)
        N_samps=length(series)
        println(f_output,"$f $avg $σ $N_samps $av $(ENV["PWD"])/histograms/histogram_$f")
        # compute histogram
        m,M = avg - n_σ*σ, avg + n_σ*σ
        dx=(M-m)/n_bins
        hist=zeros(n_bins)
        ts = (M .- series) / (M-m)
        ts = ts[0 .< ts .< 1]
        is = ceil.(Int64,ts*n_bins)
        #clamp!(is,1,n_bins)
        println(m," ", M, " ",minimum(is)," ",maximum(is))
        map(i-> hist[i]+=1,is)

        hist_file=open(joinpath("histograms","histogram_$(n_bins)_$f"),"w")
        println(hist_file,"$m $M $n_bins $avg $σ")
        print(hist_file,join(hist,", "))
        close(hist_file)
    catch e
        if isa(e,InexactError)
            println(f_output,"$f NaN NaN NaN NULL")
            println("NaN out!")
        end
    end

end

close(f_output)
#φ(m)=1/(norm*m) => φ'(m)=-1/(norm * m^2)
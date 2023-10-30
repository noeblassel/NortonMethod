using Plots,LinearAlgebra

function train_model(X,Y,N_terms;lr=1.0,grad_tol=1e-7,max_iter=100000)
    M=length(Y)
    λ=abs.(randn(N_terms))
    A=randn(N_terms)
    F=randn(N_terms)
    Ω=randn(N_terms)
    iter=0
    loss_hist=Float64[]
    ∂λ_hist=Vector{Float64}[]
    ∂A_hist=Vector{Float64}[]
    ∂F_hist=Vector{Float64}[]
    ∂Ω_hist=Vector{Float64}[]
    while iter < max_iter
        F_θ(t)=sum(A[i]*exp(-λ[i]*t)*sin(F[i]*t+Ω[i]) for i=1:N_terms)
        Y_hat=F_θ.(X)
        ΔY = Y_hat .- Y
        J= dot(ΔY,ΔY) / 2M
        println("-------------------")
        println("Iteration: ",iter," loss: ", J)
        ∂λ = [sum(-A[j] * ΔY .* X .* exp.(-λ[j] * X) .* sin.(F[j] * X .+ Ω[j]))/M for j=1:N_terms]
        ∂A = [sum(ΔY .* exp.(-λ[j]*X) .* sin.(F[j] * X .+ Ω[j]))/M for j=1:N_terms]
        ∂F = [sum(A[j] * ΔY .* X .* exp.(-λ[j] * X) .* cos.(F[j] * X .+ Ω[j]))/M for j=1:N_terms]
        ∂Ω = [sum(A[j] * ΔY .* exp.(-λ[j] * X) .* cos.(F[j] * X .+ Ω[j]))/M for j=1:N_terms]

        push!(loss_hist,J)

        push!(∂λ_hist,abs.(∂λ))
        push!(∂A_hist,abs.(∂A))
        push!(∂F_hist,abs.(∂F))
        push!(∂Ω_hist,abs.(∂Ω))

        if (all(abs.(∂λ) .<grad_tol)) && (all(abs.(∂A) .<grad_tol)) && (all(abs.(∂F) .<grad_tol)) && (all(abs.(∂Ω) .<grad_tol))
            break
        end

        λ .-= lr * ∂λ
        A .-= lr * ∂A
        F .-= lr * ∂F
        Ω .-= lr * ∂Ω
        println("-------------------")
        iter+=1
    end

    if iter == max_iter 
        println("Training has not converged")
    end

    return λ,A,F,Ω,loss_hist,∂λ_hist,∂A_hist,∂F_hist,∂Ω_hist
end


f_thevenin="acs_thevenin.txt"
f_norton="acs_norton.txt"

_,lines_thevenin...=map(split,readlines(f_thevenin))
_,lines_norton...=map(split,readlines(f_norton))

if !isdir("acs") mkdir("acs") end

t_fit_end=1.0

dt=1e-3
n_step_fit_end=round(Int64,t_fit_end/dt)

X=collect(Float64,dt*(0:n_step_fit_end-1))

N_terms=1

for (lt,ln)=zip(lines_thevenin,lines_norton)
    fnt,mt,σt,St=lt
    fnn,mn,σn,Sn=ln

    (ηt,Nxt,Nyt,Nzt)=match(r"^thevenin_response_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",fnt)
    (rn,Nxn,Nyn,Nzn)=match(r"^norton_forcing_SINUSOIDAL_(.+)_(\d+)_(\d+)_(\d+).out$",fnn)

    mt,σt,mn,σn=parse.(Float64,[mt,σt,mn,σn])
    Nxt,Nyt,Nzt,Nxn,Nyn,Nzn=parse.(Int64,[Nxt,Nyt,Nzt,Nxn,Nyn,Nzn])

    Nt=Nxt*Nyt*Nzt
    Nn=Nxn*Nyn*Nzn

    Ct=parse.(Float64,split(St,","))
    Cn=parse.(Float64,split(Sn,","))

    t_range=(1:n_step_fit_end)*dt
    plot(t_range,Ct[1:n_step_fit_end],label="thévenin",yaxis=:log)
    plot!(t_range,Cn[1:n_step_fit_end],label="norton",yaxis=:log)

    λ,A,F,Ω,loss_hist,_...=train_model(X,Ct[1:n_step_fit_end],N_terms)
    println("$λ")
    Cn_hat(t)=sum(A[i]*exp(-λ[i]*t)*sin(F[i]*t+Ω[i]) for i=1:N_terms)
    plot!(Cn_hat,label="fit")
    savefig("acs/acs_$Nn.pdf")
    plot(loss_hist,yaxis=:log,label="")
    savefig("acs/loss_$Nn.pdf")
end


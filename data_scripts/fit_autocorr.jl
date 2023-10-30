using LinearAlgebra

function train_model(X,Y,N_terms;lr=0.5,grad_tol=1e-6,max_iter=1000000)
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

X=collect(Float64,0:0.1:10)
F(x)=0.5exp(-x)*sin(2x+3)+exp(-2x)*sin(x+1)
Y=F.(X)

N_terms=2
λ,A,ν,Ω,loss_hist,∂λ_hist,∂A_hist,∂F_hist,∂Ω_hist=train_model(X,Y,N_terms)

F_θ(t)=sum(A[i]*exp(-λ[i]*t)*sin(ν[i]*t+Ω[i]) for i=1:N_terms)
using Plots
plot(X,abs.(F_θ.(X)),label="",yaxis=:log)
plot!(X,abs.(Y),label="")
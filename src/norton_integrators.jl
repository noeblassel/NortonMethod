using Random, LinearAlgebra

struct NortonSVIntegrator{TF,TG}
    dt::Float64
    η::Float64
    T::Float64
    γ::Float64
    F::TF # forcing profile
    G::TG # response profile
end
NortonSVIntegrator(dt::Float64, η::Float64, T::Float64, γ::Float64, F::Function, G::Function) = NortonSVIntegrator{typeof(F),typeof(G)}(dt, η, T, γ, F, G)

function Molly.simulate!(sys::System, sim::NortonSVIntegrator, n_steps; n_threads::Integer=Threads.nthreads(), rng=Random.GLOBAL_RNG)
    force_hist=Float64[]

    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    α = exp(-sim.γ * sim.dt)
    σ = sqrt(1 - α^2)

    accels = accelerations(sys, neighbors; n_threads=n_threads)

    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)

    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    q_y = view(coords_array, 2, :)

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)

    #compute useful dot products
    FdotG = dot(F_y, G_y)

    #initialize state on constant response manifold
    λ = (sim.η-dot(G_y,v_x))/FdotG
    v_x .+= λ * F_y
    
    run_loggers!(sys, neighbors, 0; n_threads=n_threads)
    for step_n = 1:n_steps
        #B step
        sys.velocities .+= accels * sim.dt / 2
        λ_13 = (sim.η-dot(G_y,v_x))/FdotG#analytic expression for Lagrange multiplier
        v_x .+= λ_13 * F_y 

        #A step
        sys.coords .+= sys.velocities * sim.dt
        sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

        accels .= accelerations(sys, neighbors; n_threads=n_threads)

        F_y .= sim.F.(q_y)
        G_y .= sim.G.(q_y)

        FdotG = dot(F_y, G_y)

        (isnan(FdotG)) && return force_hist #abort if system NaNs out
        λ_12 = (sim.η-dot(G_y,v_x))/FdotG #correction term to reproject momenta on manifold
        v_x .+=λ_12 * F_y
        #λ_12/Δt is a good approximation to the term ∇R_q(q_t,p_t)⋅p_t/F(q_t)⋅G(q_t) forcing the dynamics to remain on the cotangent bundle
        #B step
        sys.velocities .+= accels * sim.dt / 2
        λ_23 = (sim.η-dot(G_y,v_x))/FdotG
        v_x .+= λ_23 * F_y 
        
        #O_step
        velocities .= α*sys.velocities + σ * random_velocities(sys, sim.T; rng=rng)  #equilibrium fd solution
        λ_fd=(sim.η - dot(G_y,v_x))/FdotG #analytic expression for Lagrange multiplier
        v_x .+= λ_fd * F_y

        F_ham=(λ_13+λ_23)/sim.dt
        F_ou=sim.γ*sim.η/FdotG
        F_corr= λ_12/sim.dt 

        push!(force_hist, F_ham+F_corr+F_ou)

        neighbors = find_neighbors(sys, sys.neighbor_finder,neighbors,step_n; n_threads=n_threads)
        run_loggers!(sys, neighbors, step_n; n_threads=n_threads)
    end
    return force_hist
end

struct NortonSplitting{S,E,K,F,W,TF,TG}
    dt::S
    r::E
    T::K
    γ::F
    splitting::W
    F::TF # forcing profile
    G::TG # response profile
end

function Molly.simulate!(sys::System, sim::NortonSplitting, n_steps; n_threads::Integer=Threads.nthreads(), rng=Random.GLOBAL_RNG)
    α_eff = exp(-sim.γ * sim.dt / count('O', sim.splitting))
    σ_eff = sqrt(1 - α_eff^2) #noise on velocities, not momenta

    effective_dts=[sim.dt / count(c,sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    sys.coords .= wrap_coords.(sys.coords, (sys.boundary,))
    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, neighbors, 0; n_threads=n_threads)

    accels = accelerations(sys, neighbors; n_threads=n_threads)
    
    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)

    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    q_y = view(coords_array, 2, :)

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)

    #compute useful dot product
    FdotG = dot(F_y, G_y)

    #initialize state on constant response manifold
    λ = (sim.r-dot(G_y,v_x))/FdotG
    v_x .+= λ * F_y


    λ_hist=Float64[]

    for step_n=1:n_steps
        λ_A = λ_B = λ_O = 0.0

        for (i,op)=enumerate(sim.splitting)
            if op=='A'
                sys.coords .+= sys.velocities * effective_dts[i]
                sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

                F_y .= sim.F.(q_y)
                G_y .= sim.G.(q_y)
                FdotG = dot(F_y,G_y)
                λ = (r - dot(G_y, v_x))/ FdotG #reprojection in p onto the constant response manifold
                v_x .+= λ * F_y
                λ_A += λ

            elseif op=='B'
                ( force_computation_steps[i] ) &&  ( accels .= accelerations(sys,neighbors, n_threads=n_threads) )
                sys.velocities .+= accels * effective_dts[i]
                λ = (r- dot(G_y, v_x)) /FdotG
                v_x .+= λ * F_y
                λ_B += λ
                
            elseif op=='O'
                sys.velocities .= α_eff*sys.velocities + σ_eff * random_velocities(sys,sim.T; rng=rng)
                λ = (r- dot(G_y, v_x)) /FdotG
                v_x .+= λ * F_y

                λ_O += (1-α_eff)*sim.r /FdotG # only record bounded-variation increment, which is analytically known

            end
        end

        λ_est= (λ_A + λ_B + λ_O) / sim.dt
        push!(λ_hist, λ_est)

        run_loggers!(sys,neighbors,step_n)

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n ; n_threads=n_threads)
        end
    end
    return λ_hist
end

struct GeneralizedNortonSplitting{S,K,F,W,TF,TG,Tφ,Tξ}
    dt::S
    T::K
    γ::F
    splitting::W
    F::TF # forcing profile
    G::TG # response profile
    φ::Tφ # deterministic flow of the response
    ξ::Tξ # martingale increments of the response
end

function Molly.simulate!(sys::System, sim::GeneralizedNortonSplitting, n_steps; n_threads::Integer=Threads.nthreads(), rng=Random.GLOBAL_RNG)
    α_eff = exp(-sim.γ * sim.dt / count('O', sim.splitting))
    σ_eff = sqrt(1 - α_eff^2) #noise on velocities, not momenta

    effective_dts=[sim.dt / count(c,sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    sys.coords .= wrap_coords.(sys.coords, (sys.boundary,))
    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, neighbors, 0; n_threads=n_threads)

    accels = accelerations(sys, neighbors; n_threads=n_threads)
    dW = zero(sys.velocities)

    velocities_array = reinterpret(reshape, Float64, sys.velocities)
    coords_array = reinterpret(reshape, Float64, sys.coords)
    dW_array = reinterpret(reshape, Float64,dW)
    #views into longitudinal and transverse components

    v_x = view(velocities_array, 1, :)
    q_y = view(coords_array, 2, :)
    dW_x = view(dW_array,1,: )

    #initialize F and G vectors
    F_y = sim.F.(q_y)
    G_y = sim.G.(q_y)

    #compute useful dot product
    FdotG = dot(F_y, G_y)

    r=dot(G_y, v_x) # initial response

    λ_hist=Float64[]
    r_hist=Float64[]

    for step_n=1:n_steps
        λ_A = λ_B = λ_O = 0.0
        push!(r_hist,r)
        for (i,op)=enumerate(sim.splitting)
            if op=='A'

                ## response updates
                xi=sim.ξ(effective_dts[i]/3) #martingale increment -- note update of the response in the A steps account for 1/3 of the total evolution of the response
                r=sim.φ(r,effective_dts[i]/3) #udpate with deterministic flow

                ## tentative normal A-step

                sys.coords .+= sys.velocities * effective_dts[i]
                sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

                ## compute configuration-dependent quantities

                F_y .= sim.F.(q_y)
                G_y .= sim.G.(q_y)
                FdotG = dot(F_y,G_y)

                ## compute Lagrange multiplier
                λ_bar = (r - dot(G_y, v_x))/ FdotG # finite-variation component of the multiplier
                λ_mart = xi / FdotG #martingale component
                v_x .+= (λ_bar + λ_mart) * F_y

                λ_A += λ_bar
                ## add stochastic component to the response
                r += xi

            elseif op=='B'
                ( force_computation_steps[i] ) &&  ( accels .= accelerations(sys,neighbors, n_threads=n_threads) ) # recompute accelerations if need be
                
                ## response updates
                xi=sim.ξ(effective_dts[i]/3) #martingale increment
                r=sim.φ(r,effective_dts[i]/3) #udpate with deterministic flow

                ## tentative normal B-step
                sys.velocities .+= accels * effective_dts[i]

                ## compute Lagrange multiplier
                λ_bar = (r- dot(G_y, v_x)) /FdotG # finite-variation component of the multiplier
                λ_mart = xi / FdotG # martingale component
                v_x .+= (λ_bar + λ_mart) * F_y
                λ_B += λ_bar

            elseif op=='O'
                ## response updates
                xi=sim.ξ(effective_dts[i]/3) #martingale increment
                r=sim.φ(r,effective_dts[i]/3) #udpate with deterministic flow

                ## Tentative O-step
                dW .= σ_eff * random_velocities(sys,sim.T; rng=rng) #noise

                sys.velocities .*= α_eff #dissipation -- necessary to add fluctuation noise later to compute the finite-variation component of the multiplier

                ## compute Lagrange multiplier
                λ_bar = (r -dot(G_y, v_x))/FdotG # finite-variation component
                λ_mart = (xi - dot(G_y, dW_x)) # martingale component
                sys.velocities .+= dW
                v_x .+= (λ_bar + λ_mart) * F_y

                λ_O += λ_bar
            end
        end

        λ_est= (λ_A + λ_B +λ_O) / sim.dt

        push!(λ_hist, λ_est)
        run_loggers!(sys,neighbors,step_n)

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n ; n_threads=n_threads)
        end
    end
    
    
    return λ_hist
end

struct NortonSplittingColorDrift{S,E,K,F,W,TF}
    dt::S
    r::E
    T::K
    γ::F
    splitting::W
    F::TF # forcing profile
end

function NortonSplittingColorDrift(dt,r,T,γ,N,splitting)
    S,E,K,F,W = typeof.([dt,r,T,γ,splitting])
    cd = ones(N)
    cd[1:2:end] .= -1
    cd /= sqrt(N)

    return NortonSplittingColorDrift{S,E,K,F,W,typeof(cd)}(dt,r,T,γ,splitting,cd)
end

function Molly.simulate!(sys::System, sim::NortonSplittingColorDrift, n_steps; n_threads::Integer=Threads.nthreads(), rng=Random.GLOBAL_RNG)
    α_eff = exp(-sim.γ * sim.dt / count('O', sim.splitting))
    σ_eff = sqrt(1 - α_eff^2) #noise on velocities, not momenta

    effective_dts=[sim.dt / count(c,sim.splitting) for c in sim.splitting]

    forces_known = true
    force_computation_steps = Bool[]

    occursin(r"^.*B[^B]*A[^B]*$", sim.splitting) && (forces_known = false) #determine the need to recompute accelerations before B steps

    for op in sim.splitting
        if op == 'O'
            push!(force_computation_steps, false)
        elseif op == 'A'
            push!(force_computation_steps, false)
            forces_known = false
        elseif op == 'B'
            if forces_known
                push!(force_computation_steps, false)
            else
                push!(force_computation_steps, true)
                forces_known = true
            end
        end
    end

    sys.coords .= wrap_coords.(sys.coords, (sys.boundary,))
    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, neighbors, 0; n_threads=n_threads)

    accels = accelerations(sys, neighbors; n_threads=n_threads)
    
    velocities_array = reinterpret(reshape, Float64, sys.velocities)

    #views into longitudinal component

    v_x = view(velocities_array, 1, :)

    F = sim.F # F = G = constant

    #initialize state on constant response manifold
    λ = (sim.r-dot(F,v_x))
    v_x .+= λ * F


    λ_hist=Float64[]

    for step_n=1:n_steps
        λ_A = λ_B = λ_O = 0.0

        for (i,op)=enumerate(sim.splitting)
            if op=='A'
                sys.coords .+= sys.velocities * effective_dts[i]
                sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

                λ = (r - dot(F, v_x)) #reprojection in p onto the constant response manifold
                v_x .+= λ * F
                λ_A += λ

            elseif op=='B'
                ( force_computation_steps[i] ) &&  ( accels .= accelerations(sys,neighbors, n_threads=n_threads) )
                sys.velocities .+= accels * effective_dts[i]
                λ = (r- dot(F, v_x))
                v_x .+= λ * F
                λ_B += λ
                
            elseif op=='O'
                sys.velocities .= α_eff*sys.velocities + σ_eff * random_velocities(sys,sim.T; rng=rng)
                λ = (r- dot(F, v_x))
                v_x .+= λ * F

                λ_O += (1-α_eff)*sim.r # only record bounded-variation increment, which is analytically known

            end
        end

        λ_est= (λ_A + λ_B + λ_O) / sim.dt
        push!(λ_hist, λ_est)

        run_loggers!(sys,neighbors,step_n)

        if step_n != n_steps
            neighbors = find_neighbors(sys, sys.neighbor_finder, neighbors, step_n ; n_threads=n_threads)
        end
    end
    return λ_hist
end
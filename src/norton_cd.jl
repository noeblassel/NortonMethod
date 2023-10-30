using Molly, LinearAlgebra

include("utils.jl")
include("norton_integrators.jl")

println("Usage: T ρ dt γ r t_equilibration n_iter_sim N scheme")

T=parse(Float64,ARGS[1]) # 2.5
ρ=parse(Float64,ARGS[2]) # 1.0
dt=parse(Float64,ARGS[3]) # 1e-3
γ=parse(Float64,ARGS[4]) # 1.0
r=parse(Float64,ARGS[5]) 

t_eq=parse(Float64,ARGS[6])
n_iter_sim=parse(Int64,ARGS[7])

N=parse(Int64,ARGS[8])
splitting=ARGS[9]


L=cbrt(N/ρ)
box_size = CubicBoundary(L,L,L)
r_c = 2.0

simulator=NortonSplittingColorDrift(dt,r,T,γ,N,splitting)
simulator_eq=NortonSplittingColorDrift(1e-4,r,T,γ,N,splitting)

coords=place_atoms(N,box_size;min_dist=0.5)
atoms = [Atom(σ=2^(-1/6), ϵ=1.0, mass=1.0) for i in 1:N]
velocities = [velocity(1.0, T, 1.0) for i = 1:N]

inter=LennardJones(force_units=NoUnits,energy_units=NoUnits,cutoff=ShiftedForceCutoff(r_c),nl_only=true)
steps_nf=5
nf = (6r_c > L) ? DistanceNeighborFinder(nb_matrix=trues(N,N),n_steps=steps_nf,dist_cutoff=r_c) : CellListMapNeighborFinder(nb_matrix=trues(N,N),n_steps=steps_nf,dist_cutoff=r_c,unit_cell=box_size)

n_steps_eq=floor(Int64,t_eq/1e-4)
n_steps_sim=floor(Int64,t_eq/dt)

sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0)
simulate!(sys,simulator_eq,n_steps_eq)

sys= System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0)

for i=1:n_iter_sim
    println("iteration $i")
    force = simulate!(sys,simulator,n_steps_sim)

    if any(isnan.(force))
        throw(ErrorException("NaN values for the forcing"))
        exit(1)
    end

    f=open("thermo_results/norton_forcing_colordrift_$(N)_$(ρ)_$(T)_$(r).out","a")
    write(f,force)
    close(f)
end
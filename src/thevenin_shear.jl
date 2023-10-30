using Molly, LinearAlgebra

include("utils.jl")
println("Usage: T ρ dt γ η forcing_type=SINUSOIDAL|LINEAR|CONSTANT t_equilibration n_iter_sim Nx scheme cutoff_radius y_ratio z_ratio")

T=parse(Float64,ARGS[1])
ρ=parse(Float64,ARGS[2])
dt=parse(Float64,ARGS[3])
γ=parse(Float64,ARGS[4])
η=parse(Float64,ARGS[5])
forcing_type=ARGS[6]
t_eq=parse(Float64,ARGS[7])
n_iter_sim=parse(Int64,ARGS[8])

Nx=parse(Int64,ARGS[9])
splitting=ARGS[10]
r_c=parse(Float64,ARGS[11])
y_ratio = parse(Float64,ARGS[12])
z_ratio = parse(Float64,ARGS[13])

Ny=round(Int64,Nx*y_ratio)
Nz=round(Int64,Nx*z_ratio)

N=Nx*Ny*Nz

Lx=Nx/cbrt(ρ)
Ly=Lx*y_ratio
Lz=Lx*z_ratio

box_size = CubicBoundary(Lx,Ly,Lz)

n_steps_neighbors = 20


struct NEMD_longitudinal_forcing{F}
    forcing::F #transverse profile of the forcing
    η::Float64
end

function Molly.forces(inter::NEMD_longitudinal_forcing{F},s::System,neighbors=nothing) where {F}
    f=zero(s.velocities)
    f_x=view(reinterpret(reshape,Float64,f),1,:)
    q_y=view(reinterpret(reshape,Float64,s.coords),2,:)
    f_x .= inter.forcing.(q_y)
    return forcing.η*f
end

function fourier_response(s::System,args...;kwargs...)
    p_x=view(reinterpret(reshape,Float64,s.velocities),1,:)
    q_y=view(reinterpret(reshape,Float64,s.coords),2,:)
    Ly=s.boundary.side_lengths[2]
    N=length(s)
    return imag(dot(p_x,exp.(2im*π*q_y/Ly))/N)
end

sinus_forcing=(y-> sin(2π*y/Ly))
constant_forcing=(y -> (y<Ly/2) ? 1 : -1)
linear_forcing=(y -> (y<Ly/2) ? 4*(y-Ly/4)/Ly : 4*(3Ly/4-y)/Ly)

forcing=NEMD_longitudinal_forcing(sinus_forcing,η)

if forcing_type=="LINEAR"
    forcing=NEMD_longitudinal_forcing(linear_forcing,η)
elseif forcing_type=="CONSTANT"
    forcing=NEMD_longitudinal_forcing(constant_forcing,η)
end

nf = (3.6r_c < min(Lx,Ly,Lz)) ? CellListMapNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff= 1.2r_c,unit_cell=box_size) : DistanceNeighborFinder(nb_matrix=trues(N,N),n_steps=n_steps_neighbors,dist_cutoff=1.2r_c)

atoms=[Atom(index=i,ϵ=1.0,σ=1.0,mass=1.0) for i=1:N]
coords=place_atoms_on_3D_lattice(Nx,Ny,Nz,box_size)
velocities=[velocity(1.0,T,1.0) for i=1:N]
inter=LennardJones(cutoff=ShiftedForceCutoff(r_c),nl_only=true,force_units=NoUnits,energy_units=NoUnits)

n_steps_eq=Int64(floor(t_eq/dt))

sim=LangevinSplitting(dt=dt,friction=γ,temperature=T,splitting=splitting;remove_CM_motion=false)
sys=System(atoms=atoms,coords=coords,velocities=velocities,pairwise_inters=(inter,),general_inters=(forcing,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0)

simulate!(sys,sim,n_steps_eq)
sys=System(atoms=atoms,coords=sys.coords,velocities=sys.velocities,pairwise_inters=(inter,),general_inters=(forcing,),boundary=box_size,neighbor_finder=nf,force_units=NoUnits,energy_units=NoUnits,k=1.0,loggers=(fourier=GeneralObservableLogger(fourier_response,Float64,1),))#temp=TemperatureLogger(Float64,1)))

for i=1:n_iter_sim
    simulate!(sys,sim,n_steps_eq)
    f=open("thermo_results/thevenin_response_$(forcing_type)_$(η)_$(Nx)_$(Ny)_$(Nz).out","a")
    write(f,values(sys.loggers.fourier))
    close(f)
    empty!(sys.loggers.fourier.history)
    #= f=open("thermo_results/thevenin_temp_$(forcing_type)_$(η)_$(N).out","a")
    write(f,values(sys.loggers.temp))
    close(f)
    empty!(sys.loggers.temp.history) =#
end
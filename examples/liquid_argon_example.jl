using Molly, LinearAlgebra, StaticArrays, Random, AtomsCalculators, DelimitedFiles

include(joinpath(@__DIR__, "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "norton_integrators.jl"))

# physical parameters (LJ reduced units for argon)
# Argon: σ = 3.405 Å, ε = 0.0103 eV, m = 39.95 u
# Derived units: ε/σ ≈ 4.846 pN (force), σ/τ ≈ 1.577 Å/ps (velocity), τ ≈ 2.159 ps
# T = 85 K  → T* = kB T/ε ≈ 0.7095
# ρ = 1.41 g/cm³ → ρ* = ρ σ³/m ≈ 0.8395

T = 85 / 119.8       # reduced temperature
ρ = 0.8395            # reduced density
dt = 5e-3             # timestep (≈ 2 fs / τ ≈ 0.93e-3, using 5e-3 for stability)
γ = 1.0               # Langevin friction (= m/τ in physical units)
r_c = 2.0             # LJ cutoff radius (= 2σ, matching physical setup)
splitting = "BAOAB"

# Unit conversions: η [pN] / (ε/σ) → η*, R [Å/ps] / (σ/τ) → R*
const ε_over_σ = 4.846  # pN, LJ force unit
const σ_over_τ = 1.577  # Å/ps, LJ velocity unit

η = 5.0 / ε_over_σ     # ≈ 1.032 in reduced units (= 5 pN)
R = 0.2 / σ_over_τ     # ≈ 0.127 in reduced units (= 0.2 Å/ps)

n_steps_eq = 10_000   # equilibration steps
n_steps_prod = 100_000  # production steps

Nx = Ny = Nz = 6
N = Nx * Ny * Nz  # 216

Lx = Nx / cbrt(ρ)
Ly = Lx
Lz = Lx
box_size = CubicBoundary(Lx, Ly, Lz)

# forcing/response
F_sin(y) = sin(2π * y / Ly)
G_sin(y) = sin(2π * y / Ly) / N

atoms = [Atom(σ=1.0, ϵ=1.0, mass=1.0) for _ in 1:N]
inter = LennardJones(cutoff=ShiftedForceCutoff(r_c), use_neighbors=true)
nf = CellListMapNeighborFinder(eligible=trues(N, N), n_steps=20, dist_cutoff=1.2r_c, unit_cell=box_size)

coords_init = place_atoms_on_3D_lattice(Nx, Ny, Nz, box_size)
velocities_init = [random_velocity(1.0, T, 1.0) for _ in 1:N]

# NEMD

println("="^60)
println("NEMD simulation (η = $η)")
println("="^60)

struct NEMD_longitudinal_forcing{F}
    forcing::F
    η::Float64
end

AtomsCalculators.@generate_interface function AtomsCalculators.forces(sys, inter_nemd::NEMD_longitudinal_forcing; kwargs...)
    f = zero(sys.velocities)
    f_x = view(reinterpret(reshape, Float64, f), 1, :)
    q_y = view(reinterpret(reshape, Float64, sys.coords), 2, :)
    f_x .= inter_nemd.forcing.(q_y)
    return inter_nemd.η * f
end

function fourier_response(s::System, args...; kwargs...)
    p_x = view(reinterpret(reshape, Float64, s.velocities), 1, :)
    q_y = view(reinterpret(reshape, Float64, s.coords), 2, :)
    L = s.boundary.side_lengths[2]
    n = length(s)
    return imag(dot(p_x, exp.(2im * π * q_y / L)) / n)
end

n_bins = 20

function velocity_profile(s::System, args...; kwargs...)
    p_x = view(reinterpret(reshape, Float64, s.velocities), 1, :)
    q_y = view(reinterpret(reshape, Float64, s.coords), 2, :)
    L = s.boundary.side_lengths[2]
    profile = zeros(n_bins)
    counts = zeros(Int, n_bins)
    for i in eachindex(p_x)
        bin = clamp(floor(Int, q_y[i] / L * n_bins) + 1, 1, n_bins)
        profile[bin] += p_x[i]
        counts[bin] += 1
    end
    for b in 1:n_bins
        counts[b] > 0 && (profile[b] /= counts[b])
    end
    return profile
end

forcing_nemd = NEMD_longitudinal_forcing(F_sin, η)
sim_nemd = LangevinSplitting(dt=dt, friction=γ, temperature=T, splitting=splitting; remove_CM_motion=false)

# equilibrate

sys_nemd = System(
    atoms=atoms, coords=copy.(coords_init), velocities=copy.(velocities_init),
    pairwise_inters=(inter,), general_inters=(forcing_nemd,),
    boundary=box_size, neighbor_finder=nf,
    force_units=NoUnits, energy_units=NoUnits, k=1.0,
)

println("Equilibrating ($n_steps_eq steps)...")
simulate!(sys_nemd, sim_nemd, n_steps_eq)

# nemd production run
sys_nemd = System(
    atoms=atoms, coords=sys_nemd.coords, velocities=sys_nemd.velocities,
    pairwise_inters=(inter,), general_inters=(forcing_nemd,),
    boundary=box_size, neighbor_finder=nf,
    force_units=NoUnits, energy_units=NoUnits, k=1.0,
    loggers=(
        fourier=GeneralObservableLogger(fourier_response, Float64, 1),
        profile=GeneralObservableLogger(velocity_profile, Vector{Float64}, 1),
    ),
)

println("Production ($n_steps_prod steps)...")
simulate!(sys_nemd, sim_nemd, n_steps_prod)

R_nemd = values(sys_nemd.loggers.fourier)
mean_R = sum(R_nemd) / length(R_nemd)
println("  mean R = $mean_R (reduced) = $(mean_R * σ_over_τ) Å/ps")

profiles = values(sys_nemd.loggers.profile)
mean_profile = sum(profiles) / length(profiles)
bin_centers = [(i - 0.5) * Ly / n_bins for i in 1:n_bins]

# equilibrate norton

println("\n" * "="^60)
println("Norton simulation (R = $R)")
println("="^60)

simulator_norton = NortonSplitting(dt, R, T, γ, splitting, F_sin, G_sin)

sys_norton = System(
    atoms=atoms, coords=copy.(coords_init), velocities=copy.(velocities_init),
    pairwise_inters=(inter,), boundary=box_size, neighbor_finder=nf,
    force_units=NoUnits, energy_units=NoUnits, k=1.0,
)

println("Equilibrating ($n_steps_eq steps)...")
simulate!(sys_norton, simulator_norton, n_steps_eq)

# norton production run
sys_norton = System(
    atoms=atoms, coords=sys_norton.coords, velocities=sys_norton.velocities,
    pairwise_inters=(inter,), boundary=box_size, neighbor_finder=nf,
    force_units=NoUnits, energy_units=NoUnits, k=1.0,
)

println("Production ($n_steps_prod steps)...")
λ_norton = simulate!(sys_norton, simulator_norton, n_steps_prod)
mean_λ = sum(λ_norton) / length(λ_norton)
println("  mean λ = $mean_λ (dimensionless, expected ≈ η = $η)")

# output

outdir = joinpath(@__DIR__, "..", "results")
mkpath(outdir)

# back to physical units
const σ_Å = 3.405        # Å
const τ_ps = 2.159        # ps

R_nemd_phys = R_nemd .* σ_over_τ                      # Å/ps
bin_centers_phys = bin_centers .* σ_Å                  # Å
mean_profile_phys = mean_profile .* σ_over_τ           # Å/ps

writedlm(joinpath(outdir, "nemd_response.dlm"), R_nemd_phys)
writedlm(joinpath(outdir, "norton_lambda.dlm"), λ_norton)          # dimensionless
writedlm(joinpath(outdir, "nemd_velocity_profile.dlm"), [bin_centers_phys mean_profile_phys])

println("\nResults written to $outdir")
println("  nemd_response.dlm         — R(t) [Å/ps] from NEMD (η = $η)")
println("  norton_lambda.dlm          — λ(t) [dimensionless] from Norton (R = $(R * σ_over_τ) Å/ps)")
println("  nemd_velocity_profile.dlm  — v_x(y) [Å, Å/ps] profile ($n_bins bins)")

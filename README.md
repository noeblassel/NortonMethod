## Fixing the flux: A dual approach to computing transport coefficients
### Authors: Noé Blassel and Gabriel Stoltz

This repository contains processed data files and source code used in the production of numerical results for [this preprint](https://arxiv.org/abs/2305.08224).

### Repository structure

```
src/                        Simulation source code (Julia / Molly.jl)
  norton_integrators.jl       Custom integrators (NortonSVIntegrator, NortonSplitting,
                              GeneralizedNortonSplitting, NortonSplittingColorDrift)
  norton_shear.jl             Norton shear-viscosity simulation entry point
  norton_cd.jl                Norton color-drift simulation entry point
  thevenin_shear.jl           Thevenin (NEMD) shear-viscosity simulation entry point
  thevenin_cd.jl              Thevenin (NEMD) color-drift simulation entry point
  utils.jl                    Lattice placement helpers

sim_scripts/                Shell scripts for launching simulation campaigns
data_scripts/               Post-processing: histograms, autocorrelations, curve fitting
plot_scripts/               Plotting scripts for figures in the paper

data/                       Processed simulation output
  autocorrelations/           Autocorrelation functions
  histograms/                 Histogram data (Norton forcing / Thevenin response)
  summary_statistics/         Aggregated results (means, variances, transport coefficients)

examples/
  liquid_argon_example.jl     A simple run demonstrating NEMD/Norton duality in liquid argon using shear viscosity simulations
```

### Dependencies

Julia >= 1.9 with [Molly.jl](https://github.com/JuliaMolSim/Molly.jl) (tested on v0.23).
Install dependencies with:

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
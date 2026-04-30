# CalciumSpike.jl

CalciumSpike.jl bridges GCaMP experimental data with biophysical forward
models generated from spiking neural network simulations. It provides tools to
simulate fluorescence traces from ground-truth spike trains, deconvolve
observed ΔF/F signals, and quantitatively compare spike-inference algorithms
against known ground truth.

The package implements the physiological calcium model of
[Deneux et al. 2016](https://doi.org/10.1038/ncomms12190) (Nat Commun 7:12190)
alongside an analytical ΔF/F deconvolution pipeline and a Julia interface to
the MATLAB MLSpike toolbox, enabling direct comparison of both approaches.

## Features

- Physiological forward model (Deneux et al. 2016): calcium ODE, saturating
  indicator nonlinearity, Brownian baseline drift
- Analytical ΔF/F deconvolution assuming single-exponential calcium dynamics
- Noise-corrected population correlation analysis across trials
- MLSpike integration via MATLAB.jl (lazy-loaded; requires MATLAB ≤ 2022)
- Makie visualisation helpers (lazy-loaded)
- Multi-threaded population simulations via `Threads.@threads`

## Installation

```julia
]dev path/to/CalciumSpike.jl
```

To build the documentation locally:

```julia
julia --project=docs/ docs/make.jl
```

## Examples

Working scripts are in the `scripts/` folder:

| Script | Description |
|--------|-------------|
| `scripts/run_demo.jl` | Full MLSpike pipeline: generate → calibrate → estimate → visualise |
| `scripts/CalciumProcs/Ca_convolution.jl` | Forward model usage and parameter exploration |
| `scripts/CalciumProcs/Biophysics_test_parameters.jl` | Parameter sensitivity of the biophysical model |
| `scripts/CalciumProcs/parameter_exp.jl` | Systematic parameter sweep experiments |
| `scripts/CalciumProcs/MLSpikes.jl` | MLSpike estimation on synthetic data |
| `scripts/CalciumProcs/recordings.jl` | Loading and formatting experimental recordings |
| `scripts/CalciumProcs/retrieve_spikes.jl` | Spike retrieval and ground-truth construction |
| `scripts/noise_correction.jl` | Noise-corrected population correlation demo |

Scripts can be reworked into step-by-step tutorials on request.

## Quick start

```julia
using CalciumSpike

# Forward model with default GCaMP parameters
params = CaModel(τ=0.81, A=0.052, g=0.01, σ=0.02)

# Simulate a fluorescence trace from a Poisson spike train (10 Hz, 30 s)
spiketimes = cumsum(0.1 .* randexp(300))
trace = calcium_trace(spiketimes, 50.0, (0.0, 30.0); params)
# trace.t, trace.F, trace.c, trace.B

# Post-process: compute ΔF/F and deconvolve
pp = CaPostProcess()
dff, t = delta_f_over_f(trace.t, trace.F)
dec = calcium_postprocess(dff, t, pp)
```

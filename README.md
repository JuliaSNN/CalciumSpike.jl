# CalciumSpike.jl

Bridges GCaMP experimental data with biophysical forward models generated
from spiking neural network simulations. Simulates fluorescence traces from
ground-truth spike trains, deconvolves observed ΔF/F signals, and
quantitatively compares spike-inference algorithms against known ground
truth.

Implements the physiological calcium model of
[Deneux et al. 2016](https://doi.org/10.1038/ncomms12190) (Nat Commun
7:12190), an analytical ΔF/F deconvolution pipeline, and a Julia interface to
the MATLAB MLSpike toolbox, so both approaches can be run and compared
directly.

## Requirements

- MATLAB ≤ 2022 (for the MLSpike-backed functions; lazy-loaded via
  `Requires.jl`, so the rest of the package works without MATLAB installed)
- [MLspike](https://github.com/MLspike/spikes) (`spikes`) and
  [brick](https://github.com/thomasdeneux/brick) MATLAB toolboxes, cloned into
  one directory and added to the MATLAB path by `activate_MLSpike(path)`
  (`path` defaults to `src/lib/` of a development clone)
- Makie/CairoMakie for plotting (also lazy-loaded)

## Installation

```julia
]dev path/to/CalciumSpike.jl
```

## Capabilities

- Physiological forward model (Deneux et al. 2016): calcium ODE, saturating
  indicator nonlinearity, Brownian baseline drift (`CaModel`, `calcium_trace`)
- Analytical ΔF/F deconvolution assuming single-exponential calcium dynamics,
  with optional Gaussian smoothing (`CaPostProcess`, `calcium_postprocess`,
  `delta_f_over_f`)
- Noise-corrected population correlation analysis across trials
  (`noise_correction.jl`)
- MLSpike integration via MATLAB.jl: autocalibration, spike estimation
  (`autocalibration`, `spike_estimate`, `evaluate_MLspike`)
- Side-by-side comparison of the biophysical model against MLSpike and a
  phenomenological double-exponential kernel (`run_comparison`,
  `plot_comparison`, `ca_fr_correlation`)
- Spatial-map utilities for population-level analyses (`spatial_map.jl`)

### Randomness and reproducibility of the forward model

By default `calcium_trace` is random: the baseline drift (`η > 0`) and the
measurement noise (`σ > 0`) come from the global RNG, so each call gives a new
noise realisation (`Random.seed!` before the call makes it repeatable).

To opt in to bit-for-bit reproducibility pass an integer `seed`:

```julia
Fs, t = calcium_trace(spikes, 50.0, (0.0, 30.0); params, seed = 42)
```

The call builds a local `Random.Xoshiro(42)` and uses it for all of its noise,
so the same spikes, `params` and `seed` give identical output regardless of the
global RNG state, earlier calls, or the number of threads. For a
`Vector{Vector{Float32}}` input the neurons are processed serially from that one
stream; within one seed every call with the same spike input gives the same
noise realisation (common random numbers across perturbation conditions).
`seed` and `rng` are mutually exclusive. `biophysical_calcium(...; seed)` seeds
the calcium noise only, not the Poisson spike trains, and
`generate_synthetic_data(...; seed)` seeds the synthetic data.

Full API reference (Documenter.jl) lives in `docs/`; build it locally with
`julia --project=docs/ docs/make.jl`.

## Scripts

Everything under `scripts/` is a standalone, directly-runnable example or
analysis — none of it is loaded by the package itself.

### `scripts/CalciumProcs/` — forward model & MLSpike demos

| Script | Description |
|---|---|
| `biophysical_forward_model.jl` | Generate a synthetic calcium trace, run MLSpike autocalibration, plot the fit and drift |
| `one_spike_example.jl` | Minimal demo: spike raster + per-neuron calcium traces, saves `figures/one_spike_example.{png,svg}` |
| `noise_correction.jl` | Synthetic population data + noise-corrected vs. raw correlation heatmaps |
| `MLSpikes.jl` | MLSpike spike estimation on synthetic data, with and without autocalibration |
| `Ca_convolution.jl` | Interactive/exploratory: firing-rate × burstiness (β) sweep comparing the forward model's own deconvolution against MLSpike, saves `deconvolution_performance_norise.jld2`. Meant to be run cell-by-cell (`##`-delimited), not top-to-bottom |
| `MLSpike_deconvolution_comparison/parameter_exp.jl` | Systematic sensitivity analysis: noise vs. smoothing, saturation vs. noise, saturation vs. firing rate, using `ca_fr_correlation`/`plot_correlation_heatmaps` |
| `MLSpike_deconvolution_comparison/retrieve_spikes.jl` | Minimal `plot_comparison` demo (biophysical vs. double-exponential vs. MLSpike) |

### `scripts/GCaMP6s/` — GCaMP6s parameter fitting against Chen et al. 2013

| Script | Description |
|---|---|
| `gcamp6s_data.jl` | Empirical constants digitized from Chen et al. 2013 (Nature 499:295–300), Figure 1b (1 AP / 10 AP transients) and Figure 1f (half-decay time vs. AP count). The values were read off the published figures with an online plot-digitizing tool — no source image or CSV is kept in the repo |
| `gcamp_response_analysis.jl` | Analysis functions (`sim_traces`, `half_decay_time`, `half_time_per_spike`, `waveform_loss`) shared by the two scripts below |
| `gcamp6fit.jl` | Single-parameter-set fit visualization against the Chen et al. data, saves `gcamp6s_fit.{png,svg}` |
| `gcamp6s_optimization.jl` | Distributed Optuna search over the 6 `CaModel` fit parameters; needs a running cluster/distributed-worker setup and the Helix data root, not runnable standalone |
| `debug_sigma_leakage.jl` | Diagnostic: shows how Gaussian smoothing after deconvolution leaks activity across epoch boundaries, sweeps `σsmooth` |
| `fitting_protocol.md` | Design doc for the GCaMP6s fitting protocol (forward model, objective, free parameters) |
| `CalciumImagingPipeline_Deconvolve.m` | MATLAB: neuropil subtraction → ΔF/F → deconvolve → Gaussian filter, for suite2p-style recordings. This is the pipeline used on the real imaging data for the optimization project |

### `scripts/DeconvolutionVerification/` — MATLAB ↔ CalciumSpike.jl equivalence check

| Script | Description |
|---|---|
| `verify_matlab_pipeline.jl` | Generates a synthetic spike train, runs it through CalciumSpike's own forward model, then deconvolves the same synthetic trace two independent ways — the MATLAB pipeline above (called live via MATLAB.jl) and `calcium_postprocess` — and compares the two outputs (correlation + overlay plot). Not part of the module; a standalone check that the two pipelines agree |

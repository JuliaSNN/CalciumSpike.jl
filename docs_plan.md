# CalciumSpike.jl Documentation Plan

## Goal

Documenter.jl docs matching SpikingNeuralNetworks.jl style: `@autodocs` API, math equations, warnonly, no doctests.
Heavy description lives in `index.md` only. Per-page intros are one sentence max.

---

## Directory structure

```
docs/
├── Project.toml          # Documenter + CalciumSpike
├── make.jl
└── src/
    ├── index.md          # full package description + feature list + install
    ├── forward_model.md  # CaModel + Ca_model.jl API
    ├── postprocessing.md # CaPostProcess + post_processing.jl API
    ├── comparison.md     # model_comparison.jl API
    ├── noise_correction.md # noise_correction.jl API
    ├── mlspike.md        # MLSpike.jl MATLAB integration API
    ├── visualization.md  # plots.jl API
    └── api_reference.md  # full @autodocs dump
```

---

## Pages content outline

### `index.md` — Home (full description here)
- One-paragraph aim: bridge GCaMP experimental data with biophysical forward models from SNN-generated spike trains; compare MLSpike inference vs. analytical deconvolution against ground truth.
- Feature bullets:
  - Physiological forward model (Deneux et al. 2016, Nat Commun 7:12190)
  - Analytical ΔF/F deconvolution pipeline
  - Noise-corrected population correlation analysis across trials
  - MLSpike integration via MATLAB.jl (lazy-loaded, requires MATLAB ≤ 2022)
  - Makie visualisation (lazy-loaded)
- Install snippet (`]add` path)
- Section: **Examples** — points to `scripts/` folder:
  - `scripts/run_demo.jl` — full MLSpike pipeline demo
  - `scripts/CalciumProcs/Ca_convolution.jl` — forward model usage
  - `scripts/CalciumProcs/Biophysics_test_parameters.jl` — parameter sensitivity
  - `scripts/CalciumProcs/parameter_exp.jl` — parameter sweep experiments
  - `scripts/CalciumProcs/MLSpikes.jl` — MLSpike estimation
  - `scripts/CalciumProcs/recordings.jl` — experimental data loading
  - `scripts/CalciumProcs/retrieve_spikes.jl` — spike retrieval
  - `scripts/noise_correction.jl` — noise correction demo
  - Note: scripts can be reworked into proper tutorials on request.

### `forward_model.md` — Biophysical Forward Model
- One-sentence intro only.
- Math: calcium ODE (dual-exp), observation equation `F = B(1 + Ac/(1+gc)) + σε`, Brownian drift.
- `CaModel` parameter table (τ, A, g, F0, η, σ, τr) with units.
- `@autodocs` for: `CaModel`, `bin_spikes`, `calcium_dynamics`, `baseline_drift`, `fluorescence`, `calcium_trace`.

### `postprocessing.md` — Signal Post-Processing
- One-sentence intro.
- Math: ΔF/F definition, deconvolution formula `r(t) ≈ (dF/dt + F/τ)/A`.
- `CaPostProcess` parameter table (τ, A, σnoise, σsmooth).
- `@autodocs` for: `CaPostProcess`, `delta_f_over_f`, `gaussian_smooth`, `deconvolve_df_f`, `calcium_postprocess`.

### `comparison.md` — Model Comparison
- One-sentence intro.
- Text workflow: `spike train → [CaModel / double-exp] → ΔF/F → deconvolution → Pearson correlation`.
- `@autodocs` for: `ca_fr_correlation`, `run_comparison`, `biophysical_calcium`, `evaluate_deconvolution`, `isi_beta`.

### `noise_correction.md` — Noise Correction
- One-sentence intro: trial-averaged population activity with noise-correlation correction.
- Method note: splits trials into groups, computes split-half correlations, divides by geometric mean of within-condition noise correlations.
- `@autodocs` for: `split_mean`, `generate_synthetic_data`, `noise_correction`.

### `mlspike.md` — MLSpike Integration
- One-sentence intro.
- Prerequisites block: MATLAB ≤ 2022, Optimization Toolbox, Statistics & ML Toolbox, Signal Processing Toolbox.
- Setup: `activate_MLSpike(path)`.
- `@autodocs` for: `activate_MLSpike`, `autocalibration`, `estimation_params`, `spike_estimate`, `MLspike_estimate`, `evaluate_MLspike`.

### `visualization.md` — Visualization
- One-sentence intro (Makie lazy-loaded).
- `@autodocs` for: `plot_comparison`, `plot_correlation_heatmaps`, `plot_gCamp_and_deconvolution!`, `plot_spike_raster!`, `plot_spike_detection`.

### `api_reference.md` — API Reference
- `@autodocs` on `CalciumSpike`, split `:function` then `:type`. Mirrors SNN style.

---

## `make.jl`

```julia
using Documenter, CalciumSpike

pages = [
    "Home"              => "index.md",
    "Forward Model"     => "forward_model.md",
    "Post-Processing"   => "postprocessing.md",
    "Model Comparison"  => "comparison.md",
    "Noise Correction"  => "noise_correction.md",
    "MLSpike"           => "mlspike.md",
    "Visualization"     => "visualization.md",
    "API Reference"     => "api_reference.md",
]

makedocs(
    sitename = "CalciumSpike.jl",
    modules  = [CalciumSpike],
    warnonly = [:autodocs_block],
    format   = Documenter.HTML(),
    pages    = pages,
)

deploydocs(repo = "github.com/aquaresi/CalciumSpike.jl.git")
```

---

## Docstrings to add

`noise_correction.jl` already has docstrings. Other files need them added:

| File | Items |
|---|---|
| `Ca_model.jl` | `CaModel`, `bin_spikes`, `calcium_dynamics`, `baseline_drift`, `fluorescence`, `calcium_trace` |
| `post_processing.jl` | `CaPostProcess`, `delta_f_over_f`, `gaussian_smooth`, `deconvolve_df_f`, `calcium_postprocess` |
| `model_comparison.jl` | `ca_fr_correlation`, `run_comparison`, `biophysical_calcium`, `evaluate_deconvolution`, `isi_beta` |
| `MLSpike.jl` | `activate_MLSpike`, `autocalibration`, `estimation_params`, `spike_estimate`, `MLspike_estimate`, `evaluate_MLspike` |
| `plots.jl` | `plot_comparison`, `plot_correlation_heatmaps`, `plot_gCamp_and_deconvolution!`, `plot_spike_raster!`, `plot_spike_detection` |

Style: one-line summary + `# Arguments` / `# Returns`. Deneux et al. 2016 cited where relevant. No doctests.

---

## Not included

- No GitHub Actions / CI deploy
- No doctests
- No logic changes
- `deploydocs` placeholder URL — remove if no remote repo yet

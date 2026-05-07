# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.12.5
#     language: julia
#     name: julia-1.12
# ---

# %% [markdown]
# # Calcium Model Parameter Sensitivity Analysis
#
# Systematic quantification of how measurement noise, temporal filtering, and indicator saturation
# degrade the fidelity of calcium-imaging recordings. Implements the forward model from Deneux et al. (2016)
# to study three key parameter interactions:
#
# 1. **Noise vs. Smoothing**: Trade-off between stochastic corruption and temporal resolution
# 2. **Saturation vs. Noise**: Nonlinear distortion vs. measurement uncertainty
# 3. **Saturation vs. Firing Rate**: Activity-dependent saturation in fast-firing neurons
#
# ## Biophysical Model
#
# Following Deneux et al. (Nature Communications, 2016), normalized calcium dynamics with saturating fluorescence indicator:
#
# $$\frac{dc}{dt} = s(t) - \frac{c(t)}{\tau}$$
#
# $$\frac{dB}{dt} = \eta \, dW(t)$$
#
# $$F(t) = B(t) \left(1 + A \frac{c(t)}{1 + g \, c(t)}\right) + \sigma \, \varepsilon(t)$$
#
# where:
# - $c(t)$ is normalized intracellular calcium concentration (at rest $c=0$; after single AP $c=1$)
# - $s(t) = \sum_i \delta(t - t_i)$ is the spike train (Dirac sum)
# - $\tau$ is the calcium decay time constant (population mean: 0.81 s for OGB-1)
# - $B(t)$ is drifting baseline fluorescence (Brownian motion if $\eta > 0$)
# - $F(t)$ is measured fluorescence intensity
# - $A$ is relative fluorescence increase per spike (ΔF/F₀; OGB-1 mean: 5.2%)
# - $g$ is saturation parameter: $g = \frac{\Delta[Ca]_T}{[Ca]_0 + K_d}$ (inverse spike count at half-saturation; $g=0$ recovers linear limit)
# - $\eta$ is baseline drift amplitude (set to 0 for flat baseline)
# - $\sigma$ is white measurement-noise standard deviation
# - $\varepsilon(t)$ is unit-variance Gaussian noise
#
# **Signal fidelity** is quantified as Pearson correlation between post-processed calcium trace
# and ground-truth spike-driven firing rate (deconvolution removes indicator dynamics).

# %%
using DrWatson

using SpikingNeuralNetworks
using SNNModels
using CairoMakie
using ProgressBars

@load_units

include("CA_model.jl")

# %% [markdown]
# ## Numerical Implementation
#
# **Calcium integration:** Exact exponential decay between discrete time bins: $c[i+1] = c[i] e^{-\Delta t/\tau} + n_i$
# where $n_i$ is spike count in bin $i$ and $\Delta t = 1/f_s$ (here $f_s = 50$ Hz).
#
# **Deconvolution:** L2-regularized inverse filtering to recover firing rate from ΔF/F signals,
# separately for calcium dynamics and exponential double-exponential model (baseline for comparison).
#
# **Quality metric:** Kendall correlation (robust to outliers) between deconvolved and ground-truth rates.

# %% [markdown]
# ## Experimental Design
#
# **Standard conditions across all analyses:**
# - Baseline firing rate: 5 Hz (Poisson process)
# - Baseline saturation: g = 0.01 (modest nonlinearity)
# - Baseline noise: σ = 0 (noise-free, varied per analysis)
# - Baseline smoothing: σ_smooth = 100 ms (temporal filtering postprocessing)
# - Monte Carlo: 10 independent Poisson spike trains per condition
# - Duration: 50 seconds per realization
#
# Each analysis sweeps two parameters orthogonally while holding the rest fixed.
#
#

# %%

# Fixed parameters
n_realizations = 10
simulation_time = 50s
dt = 20ms

# Base biophysical model parameters
base_params = CaModel(
    τ  = 2s,
    A  = 0.2,
    F0 = 1.0,
    η  = 0.0,
)


# %% [markdown]
# ## Comparison: Biophysical vs. Phenomenological Models
#
# The following visualization compares three signal recovery pipelines applied to simulated calcium recordings:
#
# **1. Biophysical model (MLspike):** Deneux et al. forward model with estimated parameters (τ, A, g)
# **2. Exponential double-exponential:** Heuristic two-stage convolution (rise + decay), parameter-free
# **3. Ground truth:** Actual spike-driven firing rate (unobserved, used for benchmarking only)
#
# **Conditions across the two models:**
# - Firing rate: 5 Hz 
# - Saturation g: 0.01 
# - Measurement noise σ: 0.01 
# - Smoothing kernel: 100 ms
# - Realizations: 10 independent Poisson spike trains per condition
# - Simulation: 50 seconds → compute correlations on deconvolved signals
#
# Four-panel layout displays:
# - **Upper panels:** Raw calcium signals ($\Delta F/F_0$) for both models alongside deconvolved activity traces
# - **Time series correlation:** Cross-correlation matrices between all signal pairs (ΔF/F samples)
# - **Correlation heatmap:** Quantitative comparison (Pearson $r$) of deconvolved signals vs. ground truth
#
# This cell demonstrates how saturation (g ≠ 0) introduces systematic bias in exponential-model deconvolution,
# while proper biophysical inversion can partially recover spike timing despite nonlinearity.
#

# %%
include("Ca_convolution.jl")

params = CaModel(
    τ  = 2s,
    A  = .2,
    g  = 0.01,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.1,
)
σ = 10ms
input_rate = 5Hz
run_comparison(params, σ, input_rate)

# %%
include("Ca_convolution.jl")

params = CaModel(
    τ  = 2s,
    A  = .2,
    g  = 0.01,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.1,
)
σ = 100ms
input_rate = 5Hz
run_comparison(params, σ, input_rate)

# %% [markdown]
# ## Analysis 1: Noise vs. Smoothing Trade-off
#
# **Fixed parameters:** g = 0.01 (moderate saturation), firing rate = 15 Hz
#
# **Swept parameters:**
# - Noise σ: 20 log-spaced values from 0.01 to 1.0
# - Smoothing kernel σ_smooth: 10 linear values from 1 to 200 ms
#
# **Finding:** At fixed saturation, noise and smoothing compete. Too little smoothing preserves
# noise artifacts; too much destroys spike timing information. Optimal smoothing depends on noise level.

# %% [markdown]
# ### Results: Noise vs. Smoothing Heatmaps
#
# Three correlation matrices show how ground truth correlates with calcium model, 
# double-exponential model, and the models' cross-correlation.

# %%
# %% tags=["hide"]
input_rate = 5Hz
smooth_range = range(0.001f0, 200.f0, length=10)
σ_range = exp10.(range(log10(0.01), log10(1), length=20)) 

noise_smooth_every_g = Dict{Any, Any}()
for g in [0.0001, 0.001, 0.01]
    correlation_cross = zeros(length(σ_range), length(smooth_range))
    correlation_doubleexp = zeros(length(σ_range), length(smooth_range))
    correlation_calcium = zeros(length(σ_range), length(smooth_range))
    for (i, σ) in ProgressBar(enumerate(σ_range))
        params = @update base_params begin
            σ = σ
            g = g
        end
        Threads.@threads for j in eachindex(smooth_range)
            σ_smooth = smooth_range[j] |> Float32
            correlations = ca_fr_correlation(
                    input_rate,
                    params;
                    sim_time = simulation_time,
                    sr = 50Hz,
                    σ_smooth = σ_smooth
                )
            # Average over realizations
            correlation_calcium[i, j] = correlations.calcium
            correlation_doubleexp[i, j] = correlations.double_exp
            correlation_cross[i, j] = correlations.cross
            # @info "Completed σ = $(round(σ, digits=5)), smooth = $(round(smooth_range[j], digits=2)) ms: avg corr = $(round(correlation_calcium[i, j], digits=3))"
        end
    end
    noise_smooth_every_g[g] = (
        calcium = correlation_calcium,
        double_exp = correlation_doubleexp,
        cross = correlation_cross
    )
end



# %% [markdown]
# ### Very low non-linearity

# %%
@unpack calcium, double_exp, cross = noise_smooth_every_g[0.0001]
fig = plot_correlation_heatmaps(
    smooth_range ./ ms,
    σ_range,
    calcium,
    double_exp,
    cross,
    x_label = "Smoothing σ (ms)",
    y_label = "Noise σ",
    subtitle = "(g = 0.0001, 15Hz)",
    x_scale = Makie.identity,
    y_scale = log10
)
fig

# %%
@unpack calcium, double_exp, cross = noise_smooth_every_g[0.001]
fig = plot_correlation_heatmaps(
    smooth_range ./ ms,
    σ_range,
    calcium,
    double_exp,
    cross,
    x_label = "Smoothing σ (ms)",
    y_label = "Noise σ",
    subtitle = "(g = 0.001, 15Hz)",
    x_scale = Makie.identity,
    y_scale = log10
)
fig

# %% [markdown]
# ### Non-linearity comparable to $\sigma$

# %%
@unpack calcium, double_exp, cross = noise_smooth_every_g[0.01]
fig = plot_correlation_heatmaps(
    smooth_range ./ ms,
    σ_range,
    calcium,
    double_exp,
    cross,
    x_label = "Smoothing σ (ms)",
    y_label = "Noise σ",
    subtitle = "(g = 0.1, 15Hz)",
    x_scale = Makie.identity,
    y_scale = log10
)
fig

# %% [markdown]
# ## Analysis 2: Saturation vs. Noise Interaction
#
# **Fixed parameters:** smoothing σ_smooth = 50 ms, firing rate = 5 Hz
#
# **Swept parameters:**
# - Saturation g: 20 log-spaced values from 0.0001 to 0.1
# - Noise σ: 20 log-spaced values from 0.01 to 1.0
#
# **Finding:** Saturation and noise interact nonlinearly. At low saturation (g < 0.001), 
# correlation is primarily limited by noise. At high saturation (g > 0.01), nonlinear compression 
# dominates even at low noise, causing fundamental signal distortion that postprocessing cannot recover.

# %%
# %% tags=["hide"]
# Parameter ranges for g vs σ analysis
g_range_g_sigma = exp10.(range(log10(0.0001), log10(0.1), length=20))
σ_range_g_sigma = exp10.(range(log10(0.01), log10(1), length=20))

# Fixed parameters for this analysis
σ_smooth_fixed_g_sigma = Float32(50ms)
firing_rate_fixed_g_sigma = 5Hz

# Initialize correlation matrices: [g_idx, σ_idx]
corr_calcium_g_sigma = zeros(length(g_range_g_sigma), length(σ_range_g_sigma))
corr_doubleexp_g_sigma = zeros(length(g_range_g_sigma), length(σ_range_g_sigma))
corr_cross_g_sigma = zeros(length(g_range_g_sigma), length(σ_range_g_sigma))

# Sweep g and σ
for (i, g_val) in ProgressBar(enumerate(g_range_g_sigma))
    Threads.@threads for j in eachindex(σ_range_g_sigma)
        σ_val = σ_range_g_sigma[j]
        params = @update base_params begin
            g = g_val
            σ = σ_val
        end
        correlations = ca_fr_correlation(
            firing_rate_fixed_g_sigma,
            params;
            sim_time = simulation_time,
            sr = 50Hz,
            σ_smooth = σ_smooth_fixed_g_sigma
        )
        corr_calcium_g_sigma[i, j] = correlations.calcium
        corr_doubleexp_g_sigma[i, j] = correlations.double_exp
        corr_cross_g_sigma[i, j] = correlations.cross
        # @info "g=$(round(g_val, digits=5)), σ=$(round(σ_val, digits=5)): CA=$(round(correlations.calcium, digits=3)), DE=$(round(correlations.double_exp, digits=3)), Cross=$(round(correlations.cross, digits=3))"
    end
end

# %%

fig_g_sigma = plot_correlation_heatmaps(
    g_range_g_sigma,
    σ_range_g_sigma,
    corr_calcium_g_sigma,
    corr_doubleexp_g_sigma,
    corr_cross_g_sigma,
    x_label = "Saturation g",
    y_label = "Noise σ",
    subtitle = "(smoothing = 50ms, rate = 15Hz)"
)
fig_g_sigma

# %% [markdown]
# ## Analysis 3: Saturation × Firing Rate Activity Dependence
#
# **Fixed parameters:** smoothing σ_smooth = 50 ms, noise σ = 10^-1.5
#
# **Swept parameters:**
# - Saturation g: 20 log-spaced values from 0.01 to 1.0
# - Firing rate: 20 log-spaced values from 0.1 to 20 Hz
#
# **Finding:** Saturation effects are strongly activity-dependent. At low firing rates 
# (<1 Hz), calcium remains mostly in the linear regime regardless of g. At high firing rates 
# (>10 Hz), rapid calcium accumulation pushes even moderately saturating indicators (g=0.1) 
# deep into the nonlinear saturation regime, causing severe signal distortion.

# %%
# %% tags=["hide"]
# Parameter ranges for g vs firing rate analysis
g_range_g_rate = exp10.(range(log10(0.01), log10(1), length=20))
firing_rate_range_g_rate = exp10.(range(log10(0.1), log10(20), length=20))

# Fixed parameters for this analysis
σ_fixed_g_rate = Float32(10.0f0 ^ -1.5)
σ_smooth_fixed_g_rate = Float32(50ms)

# Initialize correlation matrices: [g_idx, rate_idx]
corr_calcium_g_rate = zeros(length(g_range_g_rate), length(firing_rate_range_g_rate))
corr_doubleexp_g_rate = zeros(length(g_range_g_rate), length(firing_rate_range_g_rate))
corr_cross_g_rate = zeros(length(g_range_g_rate), length(firing_rate_range_g_rate))

# Sweep g and firing rate
for (i, g_val) in ProgressBar(enumerate(g_range_g_rate))
    Threads.@threads for j in eachindex(firing_rate_range_g_rate)
        rate = firing_rate_range_g_rate[j] * Hz
        params = @update base_params begin
            g = g_val
            σ = σ_fixed_g_rate
        end
        correlations = ca_fr_correlation(
            rate,
            params;
            sim_time = simulation_time,
            sr = 50Hz,
            σ_smooth = σ_smooth_fixed_g_rate
        )
        corr_calcium_g_rate[i, j] = correlations.calcium
        corr_doubleexp_g_rate[i, j] = correlations.double_exp
        corr_cross_g_rate[i, j] = correlations.cross
        @info "g=$(round(g_val, digits=5)), rate=$(round(rate/Hz, digits=2))Hz: CA=$(round(correlations.calcium, digits=3)), DE=$(round(correlations.double_exp, digits=3)), Cross=$(round(correlations.cross, digits=3))"
    end
end

# %%

fig_g_rate = plot_correlation_heatmaps(
    g_range_g_rate,
    firing_rate_range_g_rate,
    corr_calcium_g_rate,
    corr_doubleexp_g_rate,
    corr_cross_g_rate,
    x_label = "Saturation g",
    y_label = "Firing Rate (Hz)",
    subtitle = "(smoothing = 50ms, σ = 10⁻¹·⁵)"
)
fig_g_rate

# %% [markdown]
# ## Summary: Three Dimensions of Signal Degradation
#
# This notebook systematically characterizes how different factors undermine calcium imaging fidelity:
#
# **Analysis 1 — Noise vs. Smoothing:** At a fixed saturation level (g=0.01), measurement noise 
# and temporal smoothing present a fundamental trade-off. Small amounts of smoothing cannot remove 
# noise without destroying spike timing; excessive smoothing blurs fast dynamics.
#
# **Analysis 2 — Saturation vs. Noise:** Saturation and noise are dependent sources of degradation. At low saturation, noise is less detrimental.
#
# **Analysis 3 — Saturation vs. Firing Rate:** The practical impact of saturation critically depends 
# on neuronal activity. Sparse neurons (<1 Hz) remain in the linear regime; fast-firing neurons (>10 Hz) 
# accumulate calcium rapidly, pushing even moderately saturating indicators into nonlinearity.
# <!-- 
# **Design Implications:**
# - Select indicators with g < 0.001 for fast-firing neurons
# - For slow neurons, focus on noise reduction rather than saturation
# - High photon detection efficiency reduces noise σ without adverse effects
# - Baseline fluorescence elevation can shift dynamics toward saturation—avoid photobleaching -->


# %% [markdown]
#

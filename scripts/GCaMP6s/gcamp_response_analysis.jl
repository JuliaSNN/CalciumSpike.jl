using DrWatson
using OptunaLoader
using CalciumSpike
using SNNModels
@load_units
using Statistics
using Distributions
import Interpolations as Itp

# Analysis functions for GCaMP fluorescence responses.
# Depends on: gcamp6s_data.jl (empirical constants), CalciumSpike, SNNModels.

## Simulation constants
const SIM_SR       = 50Hz
const SIM_INTERVAL = 0s:40s
const SIM_HEATUP   = 10s
const T_SPIKE      = 4000f0   # ms (= 4 s into simulation)

"""
Simulate ΔF/F traces for a batch of spike counts.
All spikes placed at T_SPIKE (deterministic, tight jitter).
"""
function sim_traces(spike_counts::AbstractVector{Int}, params::CaModel)
    spikes = [fill(T_SPIKE, n) |> Vector{Float32} for n in spike_counts]
    Fs, t  = calcium_trace(spikes, SIM_SR, SIM_INTERVAL; params)
    ΔFs, _ = delta_f_over_f(t, Fs; heatup_time = SIM_HEATUP)
    return ΔFs, t
end

"""
Time from peak to half-maximum in a single ΔF/F trace (seconds).
Returns 0 if peak is below threshold or half-point not found.
"""
function half_decay_time(F::AbstractVector, t)
    peak_idx = argmax(F)
    peak_val = F[peak_idx]
    peak_val < 1e-6 && return 0.0f0
    half = peak_val / 2f0
    idx  = findlast(>(half), F)
    idx === nothing && return 0.0f0
    return (t[idx] - t[peak_idx]) / 1000f0   # ms → s
end

"""
Mean half-decay time (s) for `n_spikes` action potentials.
Averages over 10 independent Gaussian spike trains (σ = 10 ms) around T_SPIKE.
"""
function half_time_per_spike(n_spikes::Int; params::CaModel, sr = SIM_SR, interval = SIM_INTERVAL)
    map(1:5) do _
        spikes    = [Float32.(rand(Distributions.Normal(T_SPIKE, 10), n_spikes))]
        gcamp, t  = calcium_trace(spikes, sr, interval; params)
        ΔF, t     = delta_f_over_f(t, gcamp; heatup_time = SIM_HEATUP)
        half_decay_time(ΔF[1], t)
    end |> mean
end

"""
Normalized waveform MSE between a simulated trace and empirical data.
Interpolates the simulation onto the empirical time grid (absolute ms).
Both traces normalized to their respective peaks before comparison.
"""
function waveform_loss(
    F         :: AbstractVector,
    t         :: AbstractVector,
    t_emp_s   :: AbstractVector,   # empirical times relative to spike (seconds)
    dF_emp    :: AbstractVector,
)
    F_itp   = Itp.scale(Itp.interpolate(collect(F), Itp.BSpline(Itp.Linear())), t)
    t_query = T_SPIKE .+ t_emp_s .* 1000f0   # relative s → absolute ms
    @assert all(t_query .>= t[begin]) && all(t_query .<= t[end]) "Empirical time points out of simulation range"
    sim_vals = F_itp.(t_query)
    peak_sim = maximum(abs.(sim_vals))
    peak_sim < 1e-8 && return 1.0f0
    return sqrt(mean((dF_emp .- sim_vals).^2)./peak_sim^2)
end

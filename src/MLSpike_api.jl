# MATLAB-backed API — implemented in ext/CalciumSpikeMATLABExt.jl, only
# usable once `using MATLAB` has also been loaded.

"""
    activate_MLSpike(path)

Add the MLSpike (`spikes`) and Brick (`brick`) MATLAB libraries to the MATLAB
path. Must be called once before any other MLSpike function.

!!! warning "MATLAB requirements"
    Requires MATLAB ≤ 2022 with the Optimization Toolbox, Statistics and
    Machine Learning Toolbox, and Signal Processing Toolbox.

# Arguments
- `path`: root directory containing `spikes/` and `brick/` subdirectories
  (defaults to `src/lib/`)

See also [`MLspike_estimate`](@ref), [`autocalibration`](@ref).
"""
function activate_MLSpike end

"""
    generate_spike_train(; ntrial=6, T=30.0, rate=1.0) -> Spiketimes

Generate bursty Poisson spike trains via the MATLAB `spk_gentrain` function.

# Arguments
- `ntrial`: number of repeated trials
- `T`: trial duration in seconds
- `rate`: mean firing rate in Hz

# Returns
- `Spiketimes`: vector of spike-time vectors (one per trial)

See also [`spk_calcium`](@ref), [`activate_MLSpike`](@ref).
"""
function generate_spike_train end

"""
    spk_calcium(spikes; params) -> (time_range, traces)

Generate synthetic calcium traces from spike trains using the MATLAB
`spk_calcium` function.

# Arguments
- `spikes::Spiketimes`: spike-time vectors to convert
- `params`: MATLAB `spk_calcium` parameter struct (see [`MLSpike_params`](@ref))

# Returns
- `(time_range, traces)`: time axis and vector of fluorescence traces

See also [`MLSpike_params`](@ref), [`spike_estimate`](@ref).
"""
function spk_calcium end

function MLSpike_params end

"""
    autocalibration(calcium; dt, amin, amax, taumin, taumax, saturation, σmin, σmax) -> NamedTuple

Estimate calcium indicator parameters from fluorescence traces using MATLAB
`spk_autocalibration`. Parameter bounds are log-uniform search ranges.

# Arguments
- `calcium`: vector of fluorescence traces
- `dt=0.02`: frame interval in seconds
- `amin, amax`: amplitude search range (default 0.04–0.1)
- `taumin, taumax`: decay time constant search range in seconds (default 0.4–1.6)
- `σmin, σmax`: noise search range (default 0.005–0.05)

# Returns
NamedTuple `(tau, a, sigma)` with the estimated decay time, amplitude, and noise.

See also [`estimation_params`](@ref), [`MLspike_estimate`](@ref).
"""
function autocalibration end

"""
    estimation_params(; dt, a, tau, saturation=0.1, sigma, drift_parameter) -> params

Construct a MATLAB `tps_mlspikes` parameter struct for spike estimation.
Logs all parameters to stdout. Use outputs of [`autocalibration`](@ref) for
`a`, `tau`, and `sigma`.

# Arguments
- `dt`: frame interval in seconds
- `a`: fluorescence amplitude per spike
- `tau`: calcium decay time constant in seconds
- `saturation=0.1`: indicator saturation parameter
- `sigma`: measurement noise standard deviation
- `drift_parameter`: baseline drift amplitude

# Returns
- MATLAB parameter struct passed to [`spike_estimate`](@ref)

See also [`autocalibration`](@ref), [`spike_estimate`](@ref).
"""
function estimation_params end

"""
    spike_estimate(calcium, par) -> NamedTuple

Infer spike trains from calcium traces using MATLAB `spk_est`.

# Arguments
- `calcium`: vector of fluorescence traces
- `par`: estimation parameter struct from [`estimation_params`](@ref)

# Returns
NamedTuple `(spikest, fit, drift)`:
- `spikest::Spiketimes`: estimated spike times
- `fit`: MLSpike fluorescence fit
- `drift`: estimated baseline drift

See also [`estimation_params`](@ref), [`MLspike_estimate`](@ref).
"""
function spike_estimate end

"""
    MLspike_estimate(calcium; dt=0.02, saturation=0.1, drift_parameter=0.01) -> NamedTuple

High-level pipeline: autocalibrate indicator parameters from `calcium[1]`, then
estimate spikes for all traces. Combines [`autocalibration`](@ref),
[`estimation_params`](@ref), and [`spike_estimate`](@ref).

# Arguments
- `calcium`: vector of fluorescence traces (first trace used for autocalibration)
- `dt=0.02`: frame interval in seconds
- `saturation=0.1`: indicator saturation
- `drift_parameter=0.01`: baseline drift amplitude

# Returns
NamedTuple `(spikest, fit, drift)` from [`spike_estimate`](@ref).

See also [`evaluate_MLspike`](@ref), [`autocalibration`](@ref).
"""
function MLspike_estimate end

"""
    evaluate_MLspike(synthetic_calcium, spikes, t; kwargs...) -> Vector{Float64}

Run [`MLspike_estimate`](@ref) on `synthetic_calcium` and return per-neuron
Pearson correlations between the estimated and ground-truth firing rates.

# Arguments
- `synthetic_calcium`: vector of fluorescence traces
- `spikes`: ground-truth spike times (vector of vectors)
- `t`: common time axis for firing rate evaluation
- `kwargs...`: forwarded to [`MLspike_estimate`](@ref)

# Returns
- Vector of per-neuron Pearson correlations

See also [`evaluate_deconvolution`](@ref), [`MLspike_estimate`](@ref).
"""
function evaluate_MLspike end

export activate_MLSpike, generate_spike_train, spk_calcium, MLSpike_params,
       spike_estimate, estimation_params, autocalibration, MLspike_estimate,
       evaluate_MLspike

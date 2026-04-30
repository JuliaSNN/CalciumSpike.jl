using MATLAB

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
function activate_MLSpike(path= joinpath(@__DIR__, "lib",))
    @info "Activating MLSpike..."
    @info "To activate MLSpike, we need to add:
        -the Optimization Toolbox
        -the Statistics and Machine Learning Toolbox
        -the Signal Processing Toolbox "
    @info "Please ensure you have these toolboxes installed and licensed in your MATLAB environment."
    @info "Adding MLSpike to MATLAB path..."
    @info "MLSpike path: $path"
    # Add MLspike and Brick libraries to MATLAB path
    # spikes :=
    # Add MLspike and Brick libraries to MATLAB path
    # spikes := https://github.com/MLspike/spikes
    # brick := https://github.com/thomasdeneux/brick
    spikes_path = joinpath(path, "spikes")
    mat"addpath($spikes_path)"
    brick_path = joinpath(path, "brick")
    mat"addpath($brick_path)"
end

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
function generate_spike_train(;ntrial::Int = 6, T = 30.0, rate = 1.0)
    # Generate calcium signal using MLSpike
    mat_spikes = nothing
    rate = Float64(rate)
    T = Float64(T)
    mat"""
    ntrial = $ntrial;
    T = $T;
    rate = $rate;
    $mat_spikes = spk_gentrain(rate, T, 'bursty', 'repeat', ntrial);
    """
    return mat_to_spiketimes(mat_spikes)
end

function spiketimes_to_mat(spikes::Spiketimes)
    spikes_array = length(spikes)
    mat_spikes = Matrix{Any}(undef, 1, spikes_array)
    for (i, trial_spikes) in enumerate(spikes)
        mat = zeros(1, length(trial_spikes))
        mat[1, :] .= trial_spikes./1000
        mat_spikes[i] = mat
    end
    return mat_spikes
end

function mat_to_spiketimes(mat_spikes)
    myvec(x::Real) = vec([x])
    myvec(x::Matrix) = vec(x)
    spiketimes = Spiketimes([myvec(mat_spikes[i]).*1000 for i in 1:size(mat_spikes, 2)])
    return spiketimes
end

function traces_to_mat(traces::Vector{Vector{Float64}})
    mat_traces = Matrix{Any}(undef, 1, length(traces))
    for (i, trace) in enumerate(traces)
        mat_traces[i] = reshape(trace, 1, :)
    end
    return mat_traces
end

function mat_to_traces(mat_traces)
    traces = [vec(mat_traces[i]) for i in 1:size(mat_traces, 2)]
    return traces
end

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
function spk_calcium(spikes; params)
    mat_spikes = nothing
    if isa(spikes, Spiketimes) 
        mat_spikes = spiketimes_to_mat(spikes)
    else
        error("Unsupported spikes format. Please provide a Spiketimes object.")
    end

    calcium = nothing
    mat"""
    $calcium = spk_calcium($mat_spikes, $params);
    """
    if length(spikes) == 1
        cc = Matrix{Any}(undef, 1, 1)
        cc[1, 1] = calcium
        calcium = cc
    end
    rr = 0:params["dt"]:(length(calcium[end]) - 1) * params["dt"]
    @assert length(calcium[1]) == length(rr)
    return rr, calcium |> mat_to_traces
end

spk_calcium(spikes, params) = spk_calcium(spikes; params=params)

function MLSpike_params(; saturation = 0.1, dt = 0.02, drift_parameter = [5, 0.015],  kwargs...)
    pcal = nothing
    mat"""
    amin = 0.04; amax = 0.1;
    taumin = 0.4; taumax = 1.6;
    sigmamin = 0.005; sigmamax = 0.05;

    a = amin * exp(rand(1) * log(amax/amin));
    tau = taumin * exp(rand(1) * log(taumax/taumin));
    sigma = sigmamin * exp(rand(1) * log(sigmamax/sigmamin));
    dt = 0.02;

    pcal = spk_calcium('par');

    pcal.dt = $dt;
    pcal.saturation = $saturation;
    pcal.drift.parameter = $drift_parameter;
    pcal.a = a;
    pcal.tau = tau;
    pcal.sigma = sigma;
    $pcal = pcal;
    """
    return pcal
    new_p = Dict{String, Any}(kwargs)
    return (;pcal..., new_p...)
end

export activate_MLSpike, generate_spike_train, spk_calcium, MLSpike_params

# ## Random physiological parameters (log-uniform)
# mat"""
# """

##
# ## Generate spikes
# spikes


# ## Generate calcium
# mat"""

# """

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
function autocalibration(calcium; dt = 0.02, amin = 0.04, amax = 0.1, taumin = 0.4, taumax = 1.6, saturation=0.1, σmin = 0.005, σmax = 0.05)
    @info "Running autocalibration with the following parameter ranges:"
    @info "a: [$amin, $amax]"
    @info "tau: [$taumin, $taumax]"
    @info "sigma: [$σmin, $σmax]"
    calcium_mat = traces_to_mat(calcium)
    tauest = nothing
    aest = nothing
    sigmaest = nothing
    mat"""
    pax = spk_autocalibration('par');
    pax.dt = $dt;
    pax.amin = $amin;
    pax.amax = $amax;
    pax.taumin = $taumin;
    pax.taumax = $taumax;
    pax.saturation = 0.1;
    pax.display = 'none';
    pax.mlspikepar.dographsummary = false;
    [$tauest $aest $sigmaest] = spk_autocalibration($calcium, pax);
    """
    a = isa(aest, Number) ? aest : NaN
    tau = isa(tauest, Number) ? tauest : NaN
    sigma = isa(sigmaest, Number) ? sigmaest : NaN
    @info "Autocalibration results: "
    @info "a=$(round(aest, digits=3))"
    @info "tau=$(round(tauest, digits=3))"
    @info "sigma=$(round(sigmaest, digits=3))"
    return (;tau, a, sigma)
end

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
function estimation_params(; dt::R, a::R, tau::R, saturation=0.1, sigma::R, drift_parameter::R) where {R <: Float64}
    param = nothing
    mat"""
    par = tps_mlspikes('par');
    par.dt = $dt;
    par.a = $a;
    par.tau = $tau;
    par.saturation = $saturation;
    par.finetune.sigma = $sigma;
    par.dographsummary = true;
    par.drift.parameter = $drift_parameter;
    $param = par;
    delete(par)
    """
    ## Pretty print the parameters
    @info "####################################################"
    @info "Estimation parameters:"
    @info "----------------------------------------------------"
    for (field, value) in pairs(param)
        if isa(value, Dict)
            field == "algo" && continue
            field == "special" && continue
            @info "$field:"
            for (subfield, subvalue) in pairs(value)
                @info "  $subfield: $subvalue"
            end
        end
        !isa(value, Number) && continue
        @info "$field: $value"
    end
    @info "----------------------------------------------------"
    for (field, value) in pairs(param)
        if !isa(value, Number)
            @info "$field is not a number, skipping..."
        end
    end
    @info "####################################################"

    return param
    # return par, (spikest=mat_to_spiketimes(spikest), fit=fit, drift=drift)
end

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
function spike_estimate(calcium, par)
    calcium_mat = traces_to_mat(calcium)
    spikest = nothing
    fit = nothing
    drift = nothing
    mat"""
    [spikest fit drift] = spk_est($calcium_mat, par);
    $spikest = spikest;
    $fit = fit;
    $drift = drift;
    """
    return (spikest = mat_to_spiketimes(spikest), fit= fit, drift= drift)
end

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
function MLspike_estimate(calcium; dt=0.02, saturation=0.1, drift_parameter=0.01)
    autocalib = autocalibration(calcium[1:1], dt=0.02)
    estimate_params = estimation_params(dt=dt, a=autocalib.a, tau=autocalib.tau, saturation=saturation, sigma=autocalib.sigma, drift_parameter=drift_parameter)
    CalciumSpike.spike_estimate(calcium, estimate_params)
end

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
function evaluate_MLspike(synthetic_calcium, spikes, t;  kwargs...)
    with_logger(ConsoleLogger(stderr, Error)) do
        estimate = CalciumSpike.MLspike_estimate(synthetic_calcium, kwargs...)
        fr_true, _ = firing_rate(spikes, t, interpolate=false, τ=100ms)
        fr_est, _ = firing_rate(estimate.spikest, t, interpolate=false, τ=100ms)
        return map(eachindex(estimate.spikest)) do i
            cor(fr_est[i, :], fr_true[i, :])
        end
    end
end

export spike_estimate, estimation_params, autocalibration, MLspike_estimate, evaluate_MLspike
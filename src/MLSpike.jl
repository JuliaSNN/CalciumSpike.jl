
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

function generate_spike_train(;ntrial::Int = 6, T::Float64 = 30.0, rate::Float64 = 1.0)
    # Generate calcium signal using MLSpike
    mat_spikes = nothing
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
        mat[1, :] .= trial_spikes
        mat_spikes[i] = mat
    end
    return mat_spikes
end

function mat_to_spiketimes(mat_spikes)
    myvec(x::Float64) = vec([x])
    myvec(x::Matrix{Float64}) = vec(x)
    spiketimes = Spiketimes([myvec(mat_spikes[i]) for i in 1:size(mat_spikes, 2)])
    return spiketimes
end

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
    return calcium
end

function MLSpike_params(; saturation = 0.1, dt = 0.02, drift_parameter = [5, 0.015],  kwargs...)
    pcal = nothing
    mat"""
    amin = 0.04; amax = 0.1;
    taumin = 0.4; taumax = 1.6;
    sigmamin = 0.005; sigmamax = 0.05;

    a = amin * exp(rand(1) * log(amax/amin))
    tau = taumin * exp(rand(1) * log(taumax/taumin))
    sigma = sigmamin * exp(rand(1) * log(sigmamax/sigmamin))
    dt = 0.02;

    pcal = spk_calcium('par');

    pcal.dt = $dt;
    pcal.saturation = $saturation;
    pcal.drift.parameter = $drift_parameter;
    pcal.a = a;
    pcal.tau = tau;
    pcal.sigma = sigma;
    $pcal = pcal;
    delete(pcal)
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

function autocalibration(calcium; dt = 0.02, amin = 0.04, amax = 0.1, taumin = 0.4, taumax = 1.6, saturation=0.1, σmin = 0.005, σmax = 0.05)
    tauest = nothing
    aest = nothing
    sigmaest = nothing
    mat"""
    pax = spk_autocalibration('par');
    pax.dt = dt;
    pax.amin = amin;
    pax.amax = amax;
    pax.taumin = taumin;
    pax.taumax = taumax;
    pax.saturation = 0.1;
    pax.display = 'none';
    pax.mlspikepar.dographsummary = false;
    [$tauest $aest $sigmaest] = spk_autocalibration($calcium, pax);
    delete(pax)
    """
    return (tau=tauest, a=aest, sigma=sigmaest)
end

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
    par.drift.parameter = .01;
    $param = par;
    delete(par)
    """
    ## Pretty print the parameters
    println("####################################################")
    println("Estimation parameters:")
    println("----------------------------------------------------")
    for (field, value) in pairs(param)
        if isa(value, Dict)
            field == "algo" && continue
            field == "special" && continue
            println("$field:")
            for (subfield, subvalue) in pairs(value)
                println("  $subfield: $subvalue")
            end
        end
        !isa(value, Number) && continue
        println("$field: $value")
    end
    println("----------------------------------------------------")
    for (field, value) in pairs(param)
        if !isa(value, Number)
            println("$field is not a number, skipping...")
        end
    end
    println("####################################################")

    return param
    # return par, (spikest=mat_to_spiketimes(spikest), fit=fit, drift=drift)
end

function spike_estimate(calcium, par)
    spikest = nothing
    fit = nothing
    drift = nothing
    mat"""
    calcium = $calcium;
    par = $par;
    [spikest fit drift] = spk_est(calcium, par);
    $spikest = spikest;
    $fit = fit;
    $drift = drift;
    """
    return (spikest = mat_to_spiketimes(spikest), fit= fit, drift= drift)
end

export spike_estimate, estimation_params, autocalibration
"""
    ca_fr_correlation(spike_rate, params_model, params_post; sim_time=40s, sr=50Hz, n_samples=10)
        -> NamedTuple

Estimate how faithfully a calcium fluorescence readout tracks firing rate by
simulating `n_samples` independent Poisson spike trains, converting each to
ΔF/F via [`calcium_trace`](@ref) and via a double-exponential gCaMP6 kernel,
then correlating each deconvolved signal with the ground-truth firing rate.

# Arguments
- `spike_rate`: target firing rate in Hz
- `params::CaModel`: calcium forward-model parameters
- `sim_time=40s`: total simulation duration
- `sr=50Hz`: calcium imaging sampling rate
- `σ_smooth=100ms`: Gaussian smoothing kernel width for deconvolved traces
- `n_samples=10`: number of independent trials to average over

# Returns
- `(calcium, double_exp, cross)`: mean Pearson correlations of the
  deconvolved calcium trace, deconvolved double-exponential trace, and their
  mutual cross-correlation, each averaged across trials

See also [`calcium_postprocess`](@ref), [`calcium_trace`](@ref).
"""
function ca_fr_correlation(spike_rate::Real,
    params_model::CaModel,
    params_post::CaPostProcess;
    sim_time=40s::Real, 
    sr::Real=50Hz,
    n_samples::Int = 10
    )

    # Generate Poisson spike train
    n1 = PoissonParameter(spike_rate)
    pop = Population(n1, N=n_samples)
    model = compose(pop, silent=true)
    monitor!(pop, :fire)
    sim!(model, sim_time)
    dt = 1/sr
    
    # Merge parameters with current g
    
    # Generate calcium trace
    interval = (0s, get_time(model))
    correlation  = tmap(1:n_samples) do s
        # Compute ground-truth firing rate
        fr, r_fr = firing_rate(pop, 0:dt:get_time(model), neurons=s)
        gcamp, r_gcamp = firing_rate(pop, 0:dt:get_time(model), kernel = gcamp6_kernel, τr = 100ms, τd = 2000ms, neurons=s)
        ca = calcium_trace(SNN.spiketimes(pop)[s], sr, interval; params_model)
        
        # Convert to ΔF/F
        signal, r_signal = delta_f_over_f(ca.t, ca.F, q=0.08)
        signal_doubleexp, r_double = delta_f_over_f(r_gcamp, gcamp(1, r_gcamp), q=0.08)
        
        # Apply smoothing and deconvolution
        dec_signal, r_dec = calcium_postprocess(signal, r_signal, params_post)
        dec_doubleexp, r_double = calcium_postprocess(signal_doubleexp, r_double, params_post)
        
        # Compute correlation with ground truth
        [cor(dec_signal, fr(1, r_dec)), cor(dec_doubleexp, fr(1, r_double)), cor(dec_signal, dec_doubleexp)]
    end 
    corr_calcium = mean([c[1] for c in correlation])
    corr_doubleexp = mean([c[2] for c in correlation])
    corr_cross = mean([c[3] for c in correlation])

    return (calcium=corr_calcium, double_exp=corr_doubleexp, cross=corr_cross)
end


"""
    run_comparison(params_model::CaModel, input_rate::Real, params_post::CaPostProcess) -> NamedTuple

Run a single-neuron simulation and compute ΔF/F traces from both the
double-exponential gCaMP6 kernel and the Deneux et al. 2016 calcium forward
model, then deconvolve both. Used internally by [`plot_comparison`](@ref).

# Arguments
- `params_model`: [`CaModel`](@ref) parameters for the forward model
- `input_rate`: Poisson firing rate in Hz
- `params_post`: [`CaPostProcess`](@ref) parameters for post-processing


# Returns
NamedTuple with fields:
- `double_exp`: `(signal, r, dec)` — ΔF/F, time axis, deconvolved trace (double-exp kernel)
- `gcamp`: `(signal, r, dec)` — ΔF/F, time axis, deconvolved trace (Deneux forward model)
- `fr`: `(signal, r)` — ground-truth firing rate and time axis

See also [`ca_fr_correlation`](@ref), [`plot_comparison`](@ref).
"""
function run_comparison(params_model::CaModel, input_rate::Real, params_post::Union{CaPostProcess, Nothing} = nothing; sim_time=40s, sr=50Hz,  β=0,)
    pop = Population(N=100,SNNModels.InhomogeneousPoissonParam(r0 = input_rate; β, τ=10ms))
    model = compose(;pop, silent=true)
    monitor!(pop, :fire)
    sim!(model, sim_time)

    dt = 1/sr
    interval = (0s, get_time(model))
    fr, r_fr = firing_rate(pop, dt:dt:get_time(model), neurons=1)
    spikes = spiketimes(model.pop)[1]

    # Compute double-exponential kernel response and Deneux forward model response
    double_exp, r = firing_rate(pop, dt:dt:get_time(model), pop_average=true, kernel = gcamp6_kernel, τr = 100ms, τd = 2s, neurons=1)
    gcamp = calcium_trace(spikes, sr, (interval[1], interval[end]), params=params_model);

    # Convert to ΔF/F
    signal_doubleexp, r_double = delta_f_over_f(r, double_exp)
    signal_mlspike, r_mlspike = delta_f_over_f(gcamp.t, gcamp.F)

    # Apply smoothing and deconvolution
    if !isnothing(params_post)
        dec_doubleexp, r_double = calcium_postprocess(signal_doubleexp, r_double, params_post)  
        dec_mlspike, r_mlspike = calcium_postprocess(signal_mlspike, r_mlspike, params_post)
    else
        dec_doubleexp, r_double = nothing, r_double
        dec_mlspike, r_mlspike = nothing, r_mlspike
    end

    return (
        double_exp=(;signal=signal_doubleexp, r=r_double, dec=dec_doubleexp, fluo=double_exp, t=r), 
        gcamp=(signal=signal_mlspike, r=r_mlspike, dec=dec_mlspike, fluo=gcamp.F, t=gcamp.t), 
        fr=(signal=fr(1, r_fr), r=r_fr),
        spikes
        )
end

function biophysical_calcium(params_model::CaModel, input_rate::Real, params_post::Union{CaPostProcess, Nothing} = nothing; sim_time=40s, sr=50Hz,  β=0, N=100)
    pop = Population(N=N,SNNModels.InhomogeneousPoissonParam(r0 = input_rate; β, τ=10ms))
    model = compose(;pop, silent=true)
    monitor!(pop, :fire)
    sim!(model, sim_time)

    dt = 1/sr
    interval = (0s, get_time(model))
    fr, r_fr = firing_rate(pop, dt:dt:get_time(model))
    spikes = spiketimes(model.pop)


    gcamps = Vector{Any}(undef, N)
    Fs = Vector{Vector{Float64}}(undef, N)
    signals= Vector{Vector{Float64}}(undef, N)
    decs = Vector{Vector{Float64}}(undef, N)
    ts = Vector{Vector{Float64}}(undef, N)
    dec_corrs = Vector{Float64}(undef, N)

    Threads.@threads     for n in 1:N
        gcamps[n] = calcium_trace(spikes[n], sr, (interval[1], interval[end]), params=params_model)
        ΔF, r_mlspike = delta_f_over_f(r_fr, gcamps[n].F)
        signals[n] = ΔF
        Fs[n] = Float64.(gcamps[n].F)
        ts[n] = gcamps[n].t
        if !isnothing(params_post)
            dec_mlspike = calcium_postprocess(ΔF, r_mlspike, params_post)
            decs[n] = dec_mlspike
            dec_corrs[n] = cor(dec_mlspike, fr(n, r_mlspike))
        end
    end

    return (
        gcamp=(signal=signals, r=r_fr, dec=decs, fluo=Fs, t=r_fr, dec_corr=dec_corrs), 
        fr=(signal=fr(1:N, r_fr), r=r_fr),
        spikes
        )
end

function isi_beta(input_rate::Real, β::Real; sim_time=40s, N=100)
    pop = Population(N=N,SNNModels.InhomogeneousPoissonParam(r0 = input_rate; β, τ=10ms))
    model = compose(;pop, silent=true)
    monitor!(pop, :fire)
    sim!(model, sim_time)
    spikes = spiketimes(model.pop)
    return SNNModels.ISI_CV(spikes) |> mean
end


function evaluate_deconvolution(spikes, deconvolution, t;  kwargs...)
    fr_true, _ = firing_rate(spikes, t, interpolate=false, τ=100ms)
    return map(eachindex(deconvolution)) do i
        cor(deconvolution[i, :], fr_true[i, :])
    end
end


export ca_fr_correlation, run_comparison,  biophysical_calcium, evaluate_deconvolution, isi_beta

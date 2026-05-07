using SpikingNeuralNetworks
using Makie, CairoMakie
using MATLAB
using DrWatson
using CalciumSpike
using Statistics
using ThreadTools
SNN.@load_units
CalciumSpike.activate_MLSpike()

##
## Test Calcium model and MLSpike deconvolution on synthetic data
params = CaModel(
    τ  = 2s,
    τr = .0s,
    A  = 0.2,
    g  = 0.1,
    F0 = 1.0,
    η  = 0.001,
    σ  = 0.1,
)
params_post = CaPostProcess(
    τ = 2s,
    A = 0.2,
    σsmooth = 100ms,
)
# input_rate = 5Hz
# logrange over rates to test deconvolution performance across a range of firing rates

fine_scale = 25
rates = exp10.(range(log10(0.1), log10(20), length=fine_scale)) 
betas = range(0, 2000, fine_scale)
corrs = zeros(length(rates), length(betas))
ML_est = zeros(length(rates), length(betas))
for n in eachindex(rates)
    @info "Testing rate $(rates[n]) Hz"
    for m in eachindex(betas)
        res = biophysical_calcium(params, rates[n]*Hz, sim_time=100s, β=betas[m], sr=50Hz, N=20, params_post)
        original = (;spikes= res.spikes, fluo = res.gcamp.fluo , t= res.gcamp.t) 
        corrs[n, m] = mean(res.gcamp.dec_corr)
        ML_est[n, m] = mean(evaluate_MLspike(original.fluo, original.spikes, original.t))
        @info "Beta: $(betas[m]), Deconvolution Correlation: $(corrs[n, m]), MLSpike Correlation: $(ML_est[n, m])"
    end
end
save("deconvolution_performance_norise.jld2", "rates", rates, "betas", betas, "corrs", corrs, "ML_est", ML_est)

##
data = DrWatson.load("deconvolution_performance_norise.jld2") |> dict2ntuple
@unpack rates, betas, corrs, ML_est = data
isi = tmap(b->isi_beta(10Hz, b), betas)
isi = betas
fig = Figure(size=(800, 500))
levels = 0:0.02:1
ax = Axis(fig[1, 1], xlabel="Firing Rate (Hz)", ylabel="ISI CV", title="Deconvolution Performance (Correlation)", xscale=log10, yticks=(betas, string.(round.(isi, digits=2))))
contourf!(ax, rates, betas, corrs, levels=levels, colormap=:inferno)
ax = Axis(fig[1, 2], xlabel="Firing Rate (Hz)", ylabel="ISI CV", title="ML Spike Performance (Correlation)", xscale=log10, yticks=(betas, string.(round.(isi, digits=2))))
contourf!(ax, rates, betas, ML_est, levels=levels, colormap=:inferno)
Colorbar(fig[1, 3], limits=(0,1), colormap=:inferno, label="Deconvolution Correlation")
fig
# save("deconvolution_performance_norise.png", fig)
##
res = biophysical_calcium(params, 1Hz, sim_time=100s, β=200, sr=50Hz, N=2, params_post)
original = (;spikes= res.spikes, fluo = res.gcamp.fluo , t= res.gcamp.t)
res.gcamp.dec[1]
CalciumSpike.plot_spike_detection(original, res)
est = evaluate_MLspike(original.fluo, original.spikes, original.t) |> mean

##


res.gcamp.dec[1] |> lines
res.gcamp.dec[1]
length(res.gcamp.r)


res.gcamp.dec
params_post = CaPostProcess(
    τ = -1s,
    A = -1,
    σsmooth = 100ms,
)
dec_signal, r_dec = calcium_postprocess(signal, r_signal, params_post)
##
calcium_gen_params = (;
    a = params.A,
    tau = params.τ,
    sigma = params.σ,
    saturation = 0.1,
    dt = 0.02,
    drift_parameter = [1, 0.015],
)
generated_spikes = CalciumSpike.generate_spike_train(ntrial=1, T=100.0, rate=input_rate/Hz)

SNNModels.ISI_CV2(generated_spikes)
SNNModels.ISI_CV2(res.spikes)

rr, calcium = CalciumSpike.spk_calcium(res.spikes,  MLSpike_params(; calcium_gen_params...))

original = (;spikes= res.spikes, fluo = calcium, t= rr)
evaluate_MLspike(original.fluo, original.spikes, original.t)


lines(betas, isi)
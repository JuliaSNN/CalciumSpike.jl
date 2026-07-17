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

res = biophysical_calcium(params, 1Hz, sim_time=30s, β=0, sr=50Hz, N=20, params_post)
original = (;spikes= res.spikes, fluo = res.gcamp.fluo , t= res.gcamp.t) 

lines(res.gcamp.signal[1])

fig = Figure(size=(400, 500))
fig, ax, plt = raster(res.spikes, markersize=5)
ax.title = "Simulated Spike Train"
ax = Axis(fig[1, 2], xlabel="Time (s)", ylabel="Fluorescence (a.u.)", title="Simulated Calcium Signal")
for n in 1:length(original.fluo)
    x = original.fluo[n] .- mean(original.fluo[n])
    lines!(ax, original.t, x .+ n, )
end
# lines!(ax, original.t, original.fluo[1], label="Simulated Calcium Signal", color=:blue)
fig
save(joinpath(@__DIR__,"..","figures","one_spike_example.png"), fig)
save(joinpath(@__DIR__,"..","figures","one_spike_example.svg"), fig)
using SpikingNeuralNetworks
using SNNModels
using Makie, CairoMakie
using Interpolations
using Statistics

@load_units

include("CA_model.jl")
# Create a sample spiketimes data structure


# fr, r = firing_rate(pop, 0:get_time(model), pop_average=true)
# fig, ax, plt = lines(r, fr, axis=(xlabel="Time (s)", ylabel="Firing Rate (Hz)", title="Firing Rate of Poisson Population"))

##


params = CaModel(
    τ  = 2s,
    A  = .2,
    g  = 0.01,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.1,
)

σ = 1ms
input_rate = 2Hz

run_comparison(params, σ, input_rate)
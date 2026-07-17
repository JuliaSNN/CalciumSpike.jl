using SpikingNeuralNetworks
using SNNModels
using CalciumSpike
using Makie, CairoMakie
using Interpolations
using Statistics

@load_units

# Create a sample spiketimes data structure

##

params = CaModel(
    τ  = 2s,
    A  = .2,
    g  = 0.01,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.1,
)

input_rate = 2Hz

plot_comparison(params, input_rate, CaPostProcess(σsmooth = 1ms))
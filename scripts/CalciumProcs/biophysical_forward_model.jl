# %% [markdown]
# # CalciumSpike biophysical forward-model example
#
# This short example shows how to generate a synthetic calcium trace, run MLSpike
# autocalibration, estimate spike times, and visualize the fit and drift.

# %% [markdown]
# ## Load libraries

# %%
using SpikingNeuralNetworks
using Makie, CairoMakie
using MATLAB
using CalciumSpike
using Statistics
@load_units

# %% [markdown]
# ## Generate synthetic calcium

# %%
params = CaModel(
    τ  = 2s,
    τr = .1s,
    A  = 0.2,
    g  = 0.1,
    F0 = 1.0,
    η  = 0.000,
    σ  = 0.1,
)

input_rate = 2Hz
sim_time = 100s
res = run_comparison(params, input_rate, sim_time=100s, sr=50Hz)
synthetic_calcium = [Float64.(res.gcamp.fluo)]

# %% [markdown]
# ## Estimate parameters with autocalibration
#
# We run MLSpike's autocalibration and pretty-print the inferred decay and noise
# parameters for quick inspection.

# %%
CalciumSpike.activate_MLSpike()
autocalib = autocalibration(synthetic_calcium, dt=0.02)

println("Autocalibration results:")
println("  a   = ", autocalib.a)
println("  tau = ", autocalib.tau)
println("  sigma = ", autocalib.sigma)

estimate_params = estimation_params(
    dt=0.02,
    a=autocalib.a,
    tau=autocalib.tau,
    saturation=0.1,
    sigma=autocalib.sigma,
    drift_parameter=0.01,
)

synt_calcium_estimate = CalciumSpike.spike_estimate(synthetic_calcium, estimate_params)

# %% [markdown]
# ## Plot spikes, calcium trace, fit, and drift

# %%
fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Neurons", title="Spike matched")
plot_spike_raster!(ax, [res.spikes], synt_calcium_estimate.spikest)
xlims!(ax, 20, sim_time / 1000)

ax = Axis(
    fig[2, 1],
    xlabel="Time (s)",
    ylabel="Fluorescence (a.u.)",
    title="MLSpike with autocalibration",
    yaxisposition=:right,
)
lines!(ax, res.gcamp.t ./ 1000, synthetic_calcium[1], label="Synthetic calcium")

ca_fit = synt_calcium_estimate.fit[1]
lines!(ax, res.gcamp.t ./ 1000, ca_fit, label="MLSpike fit", color=:red)
ca_drift = synt_calcium_estimate.drift[1]
lines!(ax, res.gcamp.t ./ 1000, ca_drift, label="Estimated drift", color=:red, linestyle=:dash)
xlims!(ax, 20, sim_time / 1000)
axislegend(ax, position=:rb)

fig


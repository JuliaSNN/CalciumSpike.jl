using CairoMakie, Makie
include(joinpath(@__DIR__, "gcamp6s_data.jl"))
include(joinpath(@__DIR__, "gcamp_response_analysis.jl"))

##

params = CaModel(
    τ  = 1.83s,
    τr = 197ms,
    A  = 0.3,
    g  = 0.05,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.01,
    c0 = 0.23,
    n  = 2.05,
)

fig = Figure(size=(800, 400))
ax  = Axis(fig[1, 1], xlabel="Number of spikes", ylabel="Half-time (s)", xscale=log10)
lines!(ax, 1:5:100, half_time_per_spike.(1:5:100; params))
scatter!(ax, DATA_X, T_HALF_EMPIRICAL, color=:red)
ylims!(ax, 0, 5)

spikes = [Float32.(rand(Distributions.Normal(T_SPIKE, 10), AP_1))]
gcamp, t = calcium_trace(spikes, SIM_SR, SIM_INTERVAL; params)
ΔF_over_F, t = delta_f_over_f(t, gcamp; heatup_time=SIM_HEATUP)
ax = Axis(fig[1, 2], xlabel="Time (s)", ylabel="Fluorescence (a.u.)")
error = waveform_loss(ΔF_over_F[1], t, T_1AP, DF_1AP)
lines!(ax, t ./ 1000, ΔF_over_F[1], color=:red, label="1 spike (error = $(round(error, sigdigits=3)))")

spikes = [Float32.(rand(Distributions.Normal(T_SPIKE, 10), AP_10))]
gcamp, t = calcium_trace(spikes, SIM_SR, SIM_INTERVAL; params)
ΔF_over_F, t = delta_f_over_f(t, gcamp; heatup_time=SIM_HEATUP)
lines!(ax, t ./ 1000, ΔF_over_F[1], color=:blue)
error = waveform_loss(ΔF_over_F[1], t, T_10AP, DF_10AP)
lines!(ax, t ./ 1000, ΔF_over_F[1], color=:blue, label="10 spikes (error = $(round(error, sigdigits=3)))")


scatter!(ax, T_1AP  .+ T_SPIKE / 1000f0, DF_1AP,  color=:red)
scatter!(ax, T_10AP .+ T_SPIKE / 1000f0, DF_10AP, color=:blue)
axislegend(ax, position=:rt)
xlims!(ax, 3, 7)
fig




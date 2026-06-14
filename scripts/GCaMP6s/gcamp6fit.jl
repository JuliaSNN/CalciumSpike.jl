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

@nature_theme
fig = Figure(size=(500, 400))
ax  = Axis(fig[1, 1][2:5,1], xlabel="Number of spikes", ylabel="Half-time (s)", xscale=log10)
scatter!(ax, DATA_X, T_HALF_EMPIRICAL, label=" Chen et al. 2013, Nature 499, 295–300", color=:black)
lines!(ax, 1:5:100, half_time_per_spike.(1:5:100; params), label="Model prediction", linewidth=5)
ylims!(ax, 0, 5)
Legend(fig[1, 1][1,1], ax, framevisible=false, labelsize=7, tellwidth = false, title="GCaMP6s response")

spikes = [Float32.(rand(Distributions.Normal(T_SPIKE, 10), AP_1))]
gcamp, t = calcium_trace(spikes, SIM_SR, SIM_INTERVAL; params)
ΔF_over_F, t = delta_f_over_f(t, gcamp)
ax = Axis(fig[1, 2][2:5,1], xlabel="Time (s)", ylabel="Fluorescence (a.u.)")
error = waveform_loss(ΔF_over_F[1], t, T_1AP, DF_1AP)
scatter!(ax, T_1AP  .+ T_SPIKE / 1000f0, DF_1AP, color=:black)
lines!(ax, t ./ 1000, ΔF_over_F[1], label="1 spike (error = $(round(error, sigdigits=3)))", linewidth=5)

spikes = [Float32.(rand(Distributions.Normal(T_SPIKE, 10), AP_10))]
gcamp, t = calcium_trace(spikes, SIM_SR, SIM_INTERVAL; params)
scatter!(ax, T_10AP .+ T_SPIKE / 1000f0, DF_10AP, color=:black)
ΔF_over_F, t = delta_f_over_f(t, gcamp)
error = waveform_loss(ΔF_over_F[1], t, T_10AP, DF_10AP)
lines!(ax, t ./ 1000, ΔF_over_F[1], label="10 spikes (error = $(round(error, sigdigits=3)))", linewidth=5)


Legend(fig[1, 2][1,1], ax, framevisible=false, labelsize=7, tellwidth = false, title="GCaMP6s response")
xlims!(ax, 3, 7)


save(joinpath(@__DIR__, "gcamp6s_fit.png"), fig)
save(joinpath(@__DIR__, "gcamp6s_fit.svg"), fig)




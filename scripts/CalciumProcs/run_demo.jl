using DrWatson
using Makie
using CalciumSpike

@load_units
CalciumSpike.activate_MLSpike()

# Create a sample spiketimes data structure
ss = generate_spike_train(ntrial=2)
##
amin = 0.04; amax = 0.1;
taumin = 0.4; taumax = 1.6;
sigmamin = 0.005; sigmamax = 0.05;
##
calcium_gen_params = (;
    a = amin * exp(rand() * log(amax/amin)),
    tau = taumin * exp(rand() * log(taumax/taumin)),
    sigma = sigmamin * exp(rand() * log(sigmamax/sigmamin)),
    saturation = 0.1,
    dt = 0.02,
    drift_parameter = [5, 0.015],
)

params = MLSpike_params(; calcium_gen_params...)
dump(params)

calcium = spk_calcium(ss; params)
##
naive_calibration = estimation_params(dt=0.02, a=0.3, tau=1.0, saturation=0.1, sigma=0.01, drift_parameter=0.01)
estimate = spike_estimate(calcium, naive_calibration)
@unpack spikest, fit, drift = estimate
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Fluorescence (a.u.)", title="MLSpike Spike without autocalibration")
r = range(0, stop=length(calcium[1])*params["dt"], length=length(calcium[1]))
lines!(ax, r, calcium[1], label="Simulated Calcium Signal")
lines!(ax, r, fit[1], label="MLSpike Fit", color=:black)
lines!(ax, r, drift[1], label="MLSpike Fit", color=:black, linestyle=:dash)
scatter!(ax, spikest[1], fill(1, length(spikest[1])), label="Estimated Spikes", color=:red)
spikest[1]
fig


##

autocalib = autocalibration(calcium, dt=0.02)
autocalib = estimation_params(dt=0.02, a=autocalib.a, tau=autocalib.tau, saturation=0.1, sigma=autocalib.sigma, drift_parameter=0.01)

estimate = spike_estimate(calcium, autocalib)
@unpack spikest, fit, drift = estimate

ax = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Fluorescence (a.u.)", title="MLSpike Spike with autocalibration")
r = range(0, stop=length(calcium[1])*params["dt"], length=length(calcium[1]))
lines!(ax, r, calcium[1], label="Simulated Calcium Signal")
lines!(ax, r, fit[1], label="MLSpike Fit", color=:black)
lines!(ax, r, drift[1], label="MLSpike Fit", color=:black, linestyle=:dash)
scatter!(ax, spikest[1], fill(1, length(spikest[1])), label="Estimated Spikes", color=:red)
spikest[1]
fig
##
println("True parameters: \na=$(round(calcium_gen_params.a, digits=3)), \ntau=$(round(calcium_gen_params.tau, digits=3)), \nsigma=$(round(calcium_gen_params.sigma, digits=3))")
using CairoMakie, Statistics, CalciumSpike

activity = generate_synthetic_data(100, 1000, 10; num_signals=0, noise_level=5.0, baseline_bias=0.0)
not_corrected = cor(mean(activity, dims=3)[:, :, 1]', mean(activity, dims=3)[:, :, 1]')
corrected = noise_correction(activity, n_group=2)

##
fig = Figure(size=(800, 400))
ax1 = Axis(fig[1, 1], title="Noise-corrected correlation")
heatmap!(ax1, corrected, colorrange=(0, 1), colormap=:inferno)
ax2 = Axis(fig[1, 2], title="Not corrected correlation")
hm = heatmap!(ax2, not_corrected, colorrange=(0, 1), colormap=:inferno)
Colorbar(fig[1, 3], hm, label="Correlation")
fig
##
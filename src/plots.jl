using Makie
using UnPack
"""
    plot_comparison(params, σ, input_rate) -> Figure

Run the double-exponential vs. BioPhysical model forward-model comparison via
`run_comparison` and produce a four-panel figure:

1. ΔF/F overlay — double-exp vs. BioPhysical model traces (Pearson r in title)
2. Deconvolved activity overlay — double-exp vs. BioPhysical model vs. ground truth
3. Correlation matrix heatmap (3×3, NaN below diagonal)

Saves the figure to `plots/comparison_doubleexp_mlspike.svg` alongside the
module source tree and returns it.

# Arguments
- `params`: forward-model parameter struct passed to `run_comparison`
- `σ`: measurement noise standard deviation
- `input_rate`: ground-truth firing rate used to drive `run_comparison`

See also [`plot_gCamp_and_deconvolution!`](@ref), [`plot_correlation_heatmaps`](@ref).
"""
function plot_comparison(params, σ, input_rate)
    res = run_comparison(params, σ, input_rate)
    Parameters.@unpack double_exp, gcamp, fr, spikes = res

    fig = Figure(size=(900,800))
    exp_fig = fig[1, 1] = GridLayout()
    Deneux_fig = fig[2, 1] = GridLayout()

    plot_gCamp_and_deconvolution!(exp_fig, double_exp.r, spikes, double_exp.signal, double_exp.dec, "DoubleExp")
    plot_gCamp_and_deconvolution!(Deneux_fig, gcamp.r, spikes, gcamp.signal, gcamp.dec, "BioPhysical model")

    @show size(double_exp.signal), size(gcamp.signal)
    c = cor(double_exp.signal, gcamp.signal)
    ax = Axis(fig[3,1][1,1], xlabel="Time (s)", ylabel="ΔF/F (%)", title="Correlation : $(round(c, digits=2))")
    plt1 = lines!(ax, double_exp.r, double_exp.signal, color=:purple, label="simple (double exponential kernel)")
    ax = Axis(fig[3,1][1,1], xlabel="Time (s)", ylabel="ΔF/F (%)", yaxisposition=:right)
    hidexdecorations!(ax)
    plt2 = lines!(ax, gcamp.r, gcamp.signal, color=:orange, label="biophysical (Deneux et al. 2016)")
    Legend(fig[3,1][1,2], [plt1, plt2], ["simple (double exponential kernel)", "biophysical (Deneux et al. 2016)"], orientation=:vertical)


    sign_cor_plt = fig[4,1] = GridLayout()
    c = cor(double_exp.dec, gcamp.dec)
    ax = Axis(sign_cor_plt[1,1], xlabel="Time (s)", ylabel="Activity (a.u.)", title="Correlation : $(round(c, digits=2))", yaxisposition=:right)
    plt1 = lines!(ax, double_exp.r ./s , double_exp.dec, color=:purple, label="simple (double exponential kernel)")
    ax = Axis(sign_cor_plt[1,1], xlabel="Time (s)", yaxisposition=:right)
    hideydecorations!(ax)
    hidexdecorations!(ax)
    plt2 = lines!(ax, gcamp.r ./s, gcamp.dec, color=:orange, label="biophysical (Deneux et al. 2016)")
    ax = Axis(sign_cor_plt[1,1], ylabel="Firing Rate (Hz)", title="Correlation : $(round(c, digits=2))", yaxisposition=:left)
    hidexdecorations!(ax)
    plt3 = lines!(ax, fr.r ./s, fr.signal, color=:black, label="ground truth")
    axislegend(ax, position=:rt)

    conditions = [double_exp.dec, gcamp.dec, fr.signal]
    @show size(conditions[1]), size(conditions[2]), size(conditions[3])
    z = zeros(3,3)
    z .= NaN
    for i in 1:3
        for j in i+1:3
            z[i,j] = cor(conditions[i], conditions[j])
        end
    end
    ax = Axis(sign_cor_plt[1,2], xticks=(1:3, ["double exp", "BioPhysical model", "ground truth"]), yticks=(1:3, ["double exp", "BioPhysical model", "ground truth"]),title="Correlation Matrix", aspect = DataAspect(), xticklabelrotation=π/4)
    plt = heatmap!(ax, z, colorrange=(-1,1), colormap=:balance)
    Colorbar(sign_cor_plt[1,3], plt, label="Correlation", width=15)
    colsize!(sign_cor_plt, 2, 100)
    fig

    fig
end

"""
    plot_correlation_heatmaps(x_data, y_data, corr_calcium, corr_doubleexp, corr_cross;
                              x_label, y_label, subtitle, 
                              x_scale=log10, y_scale=log10)

Create a 3-panel heatmap figure showing correlations between three signal processing approaches.

# Arguments
- `x_data`: Vector of x-axis values (e.g., saturation g or noise σ)
- `y_data`: Vector of y-axis values (e.g., firing rate or smoothing kernel)
- `corr_calcium`: 2D matrix of correlations (ground truth vs. calcium model)
- `corr_doubleexp`: 2D matrix of correlations (ground truth vs. double-exponential)
- `corr_cross`: 2D matrix of cross-correlations (calcium vs. double-exp)
- `x_label::String`: Label for x-axis
- `y_label::String`: Label for y-axis
- `subtitle::String`: Subtitle describing fixed parameters
- `x_scale`: Scale function for x-axis (default: log10; use nothing for linear)
- `y_scale`: Scale function for y-axis (default: log10; use nothing for linear)

# Returns
- `fig::Figure`: CairoMakie figure displaying all three heatmaps
"""
function plot_correlation_heatmaps(x_data, y_data, corr_calcium, corr_doubleexp, corr_cross;
                                   x_label::String, y_label::String, subtitle::String, 
                                   x_scale=log10, y_scale=log10)
    
    fig = Figure(size=(1200, 400))
    
    ax1 = Axis(fig[1, 1], xlabel=x_label, ylabel=y_label, 
               xscale=x_scale, yscale=y_scale,
               title="Ground Truth - Calcium Model\n$(subtitle)")
    ax2 = Axis(fig[1, 2], xlabel=x_label, ylabel=y_label, 
               xscale=x_scale, yscale=y_scale,
               title="Ground Truth - Double Exponential\n$(subtitle)")
    ax3 = Axis(fig[1, 3], xlabel=x_label, ylabel=y_label, 
               xscale=x_scale, yscale=y_scale,
               title="Double Exp - Calcium Model\n$(subtitle)")
    
    hm1 = heatmap!(ax1, x_data, y_data, colorrange=(0, 1), corr_calcium', colormap=:viridis)
    hm2 = heatmap!(ax2, x_data, y_data, colorrange=(0, 1), corr_doubleexp', colormap=:viridis)
    hm3 = heatmap!(ax3, x_data, y_data, colorrange=(0, 1), corr_cross', colormap=:viridis)
    
    Colorbar(fig[1, 4], hm1, label="Pearson correlation")
    
    return fig
end


"""
    plot_gCamp_and_deconvolution!(layout, r, spikes, signal, dec_exp, title)

Populate a Makie `GridLayout` with a three-axis overlay panel:
- Bottom axis (left): deconvolved activity `dec_exp` vs. time `r`
- Overlay axis (right): ΔF/F `signal` on a right-side y-axis
- Spike raster axis: spike times as scatter markers at y = 0.1

Axes share the same x range (`extrema(r)`); decorations are hidden on the
raster and right axes to avoid clutter. A legend is placed in column 2 of
`layout`.

# Arguments
- `layout`: target `GridLayout` (or figure position)
- `r`: time axis vector (seconds)
- `spikes`: spike-time vector (seconds)
- `signal`: ΔF/F fluorescence trace
- `dec_exp`: deconvolved activity trace
- `title::String`: axis title

# Returns
- `signal`: the fluorescence trace (passthrough)

See also [`plot_comparison`](@ref).
"""
function plot_gCamp_and_deconvolution!(layout, r, spikes, signal, dec_exp, title)
    ax = Axis(layout[1,1], xlabel="Time (s)", ylabel="Firing Rate (Hz)", title=title)
    right_ax = Axis(layout[1,1], xlabel="Time (s)", ylabel="ΔF/F (%)",  yaxisposition = :right)
    spikes_ax = Axis(layout[1,1])
    xlims!(ax, extrema(r)./s)
    xlims!(right_ax, extrema(r)./s)
    xlims!(spikes_ax, extrema(r)./s)
    ylims!(spikes_ax, (0, 1))
    ylims!(ax, (minimum(dec_exp), quantile(dec_exp, 0.999)))
    hidespines!(spikes_ax)
    hidedecorations!(spikes_ax)
    plt1 = scatter!(spikes_ax, spikes./s, fill(0.1, length(spikes)), color=:black, label="Spikes")
    plt2 = lines!(right_ax, r./s, signal, color=:red, label="ΔF/F")
    plt3 = lines!(ax, r./s, dec_exp, color=:blue, label="Deconvolved Exp")
    Legend(layout[1,2], [plt1, plt2, plt3], ["Spikes", "ΔF/F", "Deconvolved Exp"], orientation=:vertical)
    return signal
end

"""
    plot_spike_raster!(ax, spikes, estimate=nothing) -> ax

Draw spike times as arrows on `ax`. Native spikes are drawn in the default
color; estimated spikes (if provided) are drawn in dark red at `y + 0.5`.

# Arguments
- `ax`: Makie `Axis`
- `spikes`: ground-truth spike times (vector of vectors, one per neuron)
- `estimate`: optional estimated spike times in the same format

See also [`plot_spike_detection`](@ref).
"""
function plot_spike_raster!(ax, spikes, estimate=nothing)
    neurons = length(spikes)
    for n in eachindex(spikes)
        ss = spikes[n]
        ps = Point2f.(ss./1000, fill(n, length(ss)))
        vs = Vec2f.(0, fill(0.4, length(ss)))
        arrows2d!(ax, ps, vs, minshaftlength = 0,  tiplength = 0)
    end
    if estimate !== nothing
        for n in eachindex(estimate)
            isempty(estimate[n]) && continue
            ss = estimate[n]
            ps = Point2f.(ss./1000, fill(n +0.5, length(ss)))
            vs = Vec2f.(0, fill(0.4, length(ss)))
            arrows2d!(ax, ps, vs, minshaftlength = 0,  tiplength = 0, color=:darkred)
        end
    end
    ylims!(ax, (1-0.5, neurons+1.5))
    ax.ylabel = "Neuron"
    ax.xlabel = "Time (s)"
    ax.yticks = (range(1 +0.5, length(spikes).+0.5, minimum([length(spikes), 5])),
                 string.(range(1, length(spikes), minimum([length(spikes), 5]))))
    hidedecorations!(ax,label = false, ticklabels = false, ticks = false)
    return ax
end


"""
    plot_spike_detection(original, results; neuron_index=1, kwargs...) -> Figure

Produce a four-panel diagnostic figure comparing MLSpike and deconvolution
against ground truth for a single neuron:

1. Raw fluorescence with MLSpike fit and estimated drift
2. Spike raster (true vs. estimated)
3. MLSpike estimated vs. true firing rate (zoomed, 50–80 s)
4. Deconvolution vs. true firing rate (zoomed, 50–80 s)

# Arguments
- `original`: output from [`biophysical_calcium`](@ref) containing `fluo`, `spikes`, `t`
- `results`: output from [`biophysical_calcium`](@ref) containing `gcamp.dec`, `gcamp.dec_corr`
- `neuron_index=1`: which neuron to highlight in panels 1, 3, 4
- `kwargs...`: forwarded to [`MLspike_estimate`](@ref)

See also [`plot_comparison`](@ref), [`plot_spike_raster!`](@ref).
"""
function plot_spike_detection(original, results; neuron_index=1, kwargs...)
    @unpack fluo, spikes, t = original
    @unpack dec_corr, dec = results.gcamp

    estimate = CalciumSpike.MLspike_estimate(fluo, kwargs...)
    fr_true, r = firing_rate(spikes, t, interpolate=false, τ=100ms)
    fr_est, _ = firing_rate(estimate.spikest, t, interpolate=false, τ=100ms)
    corr =  map(eachindex(estimate.spikest)) do i
        cor(fr_est[i, :], fr_true[i, :])
    end |> mean

    sim_time = maximum(t)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Fluorescence (a.u.)", title="MLSpike Spike with autocalibration", yaxisposition=:right)
    lines!(ax, original.t ./1000, original.fluo[neuron_index], label="Simulated Calcium Signal")

    ca_fit = estimate.fit[1]
    lines!(ax, original.t./1000, ca_fit, label="Simulated Calcium Signal", color=:red)
    ca_drift = estimate.drift[1]
    lines!(ax, original.t./1000, ca_drift, label="Simulated Calcium Signal", color=:red, linestyle=:dash)
    xlims!(ax, 20, sim_time/1000)
    # ax = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Fluorescence (a.u.)", title="MLSpike Spike with autocalibration")
    # # lines!(ax, res.gcamp.r./1000, res.gcamp.signal, label="Experimental Data", color=:black)
    xlims!(ax, 20, sim_time/1000)

    ax = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Neurons", title="Spike matched")
    plot_spike_raster!(ax, original.spikes, estimate.spikest)
    xlims!(ax, 20, sim_time/1000)


    ax = Axis(fig[3, 1], xlabel="Time (s)", ylabel="Firing Rate (Hz)", title="MLSpike Corr: $(round(corr, digits=2))")
    lines!(ax, r./s, fr_true[neuron_index, :], label="True firing rate", color=:black)
    lines!(ax,  r./s, fr_est[neuron_index, :], label="Estimated firing rate", color=:orange)
    xlims!(ax, 50, 80)

    ax = Axis(fig[4, 1], xlabel="Time (s)", ylabel="Firing Rate (Hz)", title="Deconvolution Corr: $(round(mean(dec_corr), digits=2))")
    lines!(ax, r./s, fr_true[neuron_index, :], label="True firing rate", color=:black)
    xlims!(ax, 50, 80)
    ax = Axis(fig[4, 1])
    lines!(ax,  r./s, dec[neuron_index], label="Estimated firing rate", color=:orange)
    xlims!(ax, 50, 80)
    axislegend(ax)
    fig
end

export plot_comparison, plot_correlation_heatmaps, plot_gCamp_and_deconvolution!, plot_spike_raster!, plot_spike_detection
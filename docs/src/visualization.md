# Visualization

Makie-based plotting helpers for fluorescence traces and spike inference results.

```@meta
CurrentModule = CalciumSpike
```

!!! warning "Lazy loading"
    Plotting functions are loaded only when `Makie.jl` (or a Makie backend) is
    present. Load Makie before CalciumSpike, or trigger loading with
    `import Makie` after `using CalciumSpike`.

## Functions

```@autodocs
Modules = [CalciumSpike]
Order   = [:function]
Filter  = t -> nameof(t) in (:plot_comparison, :plot_correlation_heatmaps,
                              Symbol("plot_gCamp_and_deconvolution!"),
                              Symbol("plot_spike_raster!"),
                              :plot_spike_detection)
```

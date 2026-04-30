"""
    CalciumSpike

Simulate calcium fluorescence traces from spike trains using the Deneux et al.
2016 physiological forward model, and interface with the MATLAB MLspike toolbox
for spike inference from observed fluorescence.

Core exports: [`CaModel`](@ref), [`calcium_trace`](@ref) and lower-level
primitives [`bin_spikes`](@ref), [`calcium_dynamics`](@ref),
[`baseline_drift`](@ref), [`fluorescence`](@ref).

MATLAB-backed functions (`activate_MLSpike`, `spike_estimate`, …) are loaded
lazily via `Requires` when MATLAB.jl is present. Makie plotting helpers are
similarly loaded on demand.
"""
module CalciumSpike

    using Logging
    using Requires
    using DrWatson

    using Random
    using Statistics
    using ThreadTools
    import Interpolations as Itp
    using Parameters
    using UnPack

    using SNNModels
    @load_units

    function __init__()
        @require MATLAB = "10e44e05-a98a-55b3-a45b-ba969058deb6" include("MLSpike.jl")
        @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("plots.jl")
    end
    include("Ca_model.jl")
    include("post_processing.jl")
    include("model_comparison.jl")

end # module CalciumSpike

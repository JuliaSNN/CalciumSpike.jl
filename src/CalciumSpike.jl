"""
    CalciumSpike

Simulate calcium fluorescence traces from spike trains using the Deneux et al.
2016 physiological forward model, and interface with the MATLAB MLspike toolbox
for spike inference from observed fluorescence.

Core exports: [`CaModel`](@ref), [`calcium_trace`](@ref) and lower-level
primitives [`bin_spikes`](@ref), [`calcium_dynamics`](@ref),
[`baseline_drift`](@ref), [`fluorescence`](@ref).

MATLAB-backed functions (`activate_MLSpike`, `spike_estimate`, …) are
implemented in the `CalciumSpikeMATLABExt` package extension, loaded
automatically when MATLAB.jl is also loaded. Makie plotting helpers are
loaded on demand via `Requires`.
"""
module CalciumSpike

    using Requires

    using Random
    using Statistics
    using ThreadTools
    using LoopVectorization
    import Interpolations as Itp
    using Parameters
    using UnPack

    using SNNModels
    @load_units

    function __init__()
        @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("plots.jl")
    end
    include("Ca_model.jl")
    include("post_processing.jl")
    include("model_comparison.jl")
    include("noise_correction.jl")
    include("spatial_map.jl")
    include("MLSpike_api.jl")

end # module CalciumSpike

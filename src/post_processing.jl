
"""
    CaPostProcess

Parameters for the full ΔF/F → deconvolution → smoothing pipeline.

# Fields
- `τ::Float32 = 2.0s`: calcium decay time constant for deconvolution
- `A::Float32 = 1.0`: ΔF/F amplitude per spike; divides `(dF/dt + F/τ)` in [`deconvolve_df_f`](@ref). Set to 1 to match Sophie/Deneux convention.
- `σnoise::Float32 = 0.1`: expected measurement noise standard deviation
- `σsmooth::Float32 = 100ms`: Gaussian smoothing kernel width applied after deconvolution
- `skewed::Symbol = :none`: smoothing kernel shape — `:none` (symmetric), `:left` (causal, past only), `:right` (anti-causal)
- `baseline_window::Float32 = 3000ms`: running-median window for F0 estimation in [`delta_f_over_f`](@ref)

See also [`calcium_postprocess`](@ref), [`deconvolve_df_f`](@ref), [`delta_f_over_f`](@ref).
"""
CaPostProcess

@kwdef struct CaPostProcess
    τ::Float32              = 2.0s
    A::Float32              = 1.0
    σnoise::Float32         = 0.1
    σsmooth::Float32        = 100ms
    skewed::Symbol          = :none
    baseline_window::Float32 = 3000f0
end

"""
    delta_f_over_f(t, F; baseline_window=3000f0, q=0.20f0) -> (ΔF/F, t)

Compute ΔF/F from a raw fluorescence trace using a sliding-window baseline.

Applies a causal running median over `baseline_window` ms, then takes the
`q`-th quantile of that smoothed trajectory as a single scalar F0.
This matches Sophie's MATLAB pipeline (fn_filt → prctile).

# Arguments
- `t`: time axis in ms
- `F`: raw fluorescence vector, same length as `t`
- `baseline_window`: running-median window in ms (default 3000 ms = 3 s)
- `q`: quantile for F0 estimate (default 0.20 = 20th percentile)

# Returns
- `(ΔF/F, t)`: `Float32` ΔF/F and the original time axis

See also [`calcium_postprocess`](@ref).
"""
function delta_f_over_f(t::T, F::Vector{R}; baseline_window=3000f0, q=0.20f0) where {T<:AbstractVector, R<:Real}
    dt = t[2] - t[1]
    w  = max(1, round(Int, baseline_window / dt))
    baseline = [median(@view F[max(1, i - w + 1):i]) for i in eachindex(F)]
    F0 = quantile(Float64.(baseline), Float64(q))
    return Float32.(@. (F - F0) / F0), t
end

function delta_f_over_f(t::T, Fs::Vector{Vector{R}}; baseline_window=3000f0, q=0.20f0) where {T<:AbstractVector, R<:Real}
    tmap(x -> delta_f_over_f(t, x; baseline_window, q)[1], Fs), t
end



"""
    deconvolve_df_f(F, frame_rate, τ, A=0.2) -> Vector

Analytically deconvolve a ΔF/F trace to recover an estimated firing rate,
assuming single-exponential calcium dynamics with decay `τ` and amplitude `A`:

    r(t) ≈ (dF/dt + F/τ) / A

The time derivative is computed as a forward difference; the last sample is
padded with zero. No smoothing is applied — call [`gaussian_smooth`](@ref)
afterwards if needed.

# Arguments
- `F`: ΔF/F trace
- `frame_rate`: imaging frame rate in Hz
- `τ`: calcium decay time constant (seconds)
- `A`: ΔF/F per spike (default 0.2)

See also [`calcium_postprocess`](@ref), [`gaussian_smooth`](@ref).
"""
function deconvolve_df_f(F::T, frame_rate::R, τ::R, A::R=0.2) where {T<:AbstractVector, R<:Real}
      dt = 1.0 / frame_rate                                            
      dF = vcat(diff(F), 0.0) 
      r = @. (dF / dt + F / τ) / A
      return r                                               
end             


"""
    calcium_postprocess(signal, r, params) -> (dec, r)

Deconvolve and smooth a ΔF/F trace in one step:
1. [`deconvolve_df_f`](@ref) with decay `τ` and amplitude `A`
2. [`gaussian_smooth`](@ref) with kernel width `σ`

The time axis `r` is regularized to a uniform step before processing.

# Arguments
- `signal`: ΔF/F trace

# Keyword Arguments
- `τ`: calcium decay time constant (seconds)
- `r`: corresponding time axis
- `σ`: Gaussian smoothing width (same units as `r`)
- `A`: ΔF/F per spike, passed to `deconvolve_df_f`

# Returns
- `(dec, r)`: deconvolved and smoothed activity, regularized time axis

See also [`deconvolve_df_f`](@ref), [`gaussian_smooth`](@ref).
"""
function calcium_postprocess(fluo::Vector{R}, r::T, params::CaPostProcess) where {T<:AbstractVector,R<:Real}
    δT = r[2] - r[1]
    r = r[1]:δT:r[end]
    dec_exp = deconvolve_df_f(fluo, 1/δT, params.τ, params.A)
    dec_exp = gaussian_smooth(r, dec_exp, params.σsmooth; skewed=params.skewed)
    return dec_exp
end

function calcium_postprocess(fluos::Vector{Vector{R}}, r::T, params::CaPostProcess) where {T<:AbstractVector,R<:Real} 
    tmap(fluo->calcium_postprocess(fluo, r, params), fluos)
end



export calcium_postprocess, gaussian_smooth, deconvolve_df_f, delta_f_over_f, CaPostProcess
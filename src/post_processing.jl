
"""
    CaPostProcess

Parameters for the ΔF/F deconvolution and smoothing pipeline.

# Fields
- `τ::Float32 = 2.0s`: calcium decay time constant used for deconvolution
- `A::Float32 = 0.2`: ΔF/F amplitude per spike passed to [`deconvolve_df_f`](@ref)
- `σnoise::Float32 = 0.1`: expected measurement noise standard deviation
- `σsmooth::Float32 = 100ms`: Gaussian smoothing kernel width applied after deconvolution

See also [`calcium_postprocess`](@ref), [`deconvolve_df_f`](@ref).
"""
CaPostProcess

@kwdef struct CaPostProcess
    τ::Float32  = 2.0s
    A::Float32 = 0.2
    σnoise::Float32 = 0.1
    σsmooth::Float32 = 100ms
end

"""
    delta_f_over_f(t, F; q=0.08) -> (ΔF/F, t)

Compute ΔF/F (%) from a raw fluorescence trace, discarding the first 20 s
as warmup. Baseline `F0` is the `q`-th quantile of the remaining signal.

# Arguments
- `t`: time axis with physical units (seconds)
- `F`: raw fluorescence vector, same length as `t`
- `q::Real=0.08`: quantile used to estimate baseline

# Returns
- `(ΔF/F, t_trimmed)`: percent ΔF/F as `Float32` and the trimmed time axis

See also [`calcium_postprocess`](@ref).
"""
function delta_f_over_f(t::AbstractVector, F::AbstractVector; heatup_time=10s)
    heatup = findall(t .> heatup_time)
    F0 = mean(F[heatup])
    return Float32.(@. (F - F0) / F0), t
end


"""
    gaussian_smooth(xs, x, sigma) -> Vector

Apply a normalized Gaussian kernel to signal `x` sampled on grid `xs`.
Kernel half-width is `3σ` (truncated); boundary bins are renormalized by
accumulated kernel weight so edge values are not biased toward zero.

# Arguments
- `xs`: sample-position grid (used only for its step size `xs[2]-xs[1]`)
- `x`: signal to smooth, length `n`
- `sigma`: Gaussian standard deviation in the same units as `xs`

# Returns
- smoothed signal, same length as `x`

See also [`calcium_postprocess`](@ref), [`deconvolve_df_f`](@ref).
"""
function gaussian_smooth(xs::RT, x::T, σ::R) where {T<:AbstractVector, R<:Real, RT<:AbstractVector}
    step_x = xs[2] - xs[1]
    half = ceil(Int, 3σ/step_x) ## 3σ cutoff for the kernel                                                          
    xs = -half:half                                                                         
    kernel = exp.(-xs.^2 ./ (2σ^2))                                                     
    kernel ./= sum(kernel)                                                                  
    n = length(x) ## length of the input signal
    out = similar(x) 
    s = 0.0
    w = 0.0
    @inbounds for i in 1:n
        s = 0.0
        w = 0.0                                                                             
        @fastmath for (j, k) in enumerate(kernel)
            idx = i + (j - half - 1)
            idx < 1 && continue
            idx > n && continue
            s += k * x[idx]
            w += k                                                                      
        end                                                                                 
        out[i] = s / w
    end
    @assert length(out) == length(x)
    return out
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
function calcium_postprocess(signal::RT, r::T, params::CaPostProcess) where {T<:AbstractVector,RT<:AbstractVector}
    δT = r[2] - r[1]
    r = r[1]:δT:r[end]
    dec_exp = deconvolve_df_f(signal, 1/δT, params.τ, params.A)
    dec_exp = gaussian_smooth(r, dec_exp, params.σsmooth)
    return dec_exp
end

export calcium_postprocess, gaussian_smooth, deconvolve_df_f, delta_f_over_f, CaPostProcess
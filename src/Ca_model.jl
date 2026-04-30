"""
    Ca_model

Forward model of calcium fluorescence from spike trains, following the
physiological model used by MLspike.

Reference:
    Deneux, T., Kaszas, A., Szalay, G., Katona, G., Lakner, T., Grinvald, A.,
    Rózsa, B., & Vanzetta, I. (2016). Accurate spike estimation from noisy
    calcium signals for ultrafast three-dimensional imaging of large
    neuronal populations in vivo. *Nature Communications*, 7, 12190.
    doi:10.1038/ncomms12190

Model equations (Deneux et al. 2016, eq. 8):

    dc/dt = s(t) - c(t) / τ
    dB/dt = η dW(t)
    F(t)  = B(t) * (1 + A * c(t) / (1 + g * c(t))) + σ ε(t)

where `c` is normalized intracellular calcium (c = 0 at rest, c = 1 after a
single AP), `B` is the drifting baseline (Brownian), and `F` the measured
fluorescence. The linear convolution limit is recovered for `g = 0`, `η = 0`.
"""


"""
    Ca_params

Default parameters for the Deneux et al. 2016 calcium forward model, set
to the population mean reported for the synthetic indicator OGB-1 in rat
barrel cortex in vivo (Deneux et al. 2016, autocalibration section, n=24
neurons): `τ = 0.81 ± 0.40 s`, `A = 5.2 ± 1.6 %`.

The saturation `g`, drift `η` and measurement noise `σ` are not given as
population values in the paper — they are either fit per cell by
autocalibration, or set to zero for linear-limit simulations. Zero is used
here as the neutral default.

Fields:
- `τ`  : calcium decay time constant (seconds). `τ = (1 + k_S + k_B) / g_e`.
         Paper OGB mean: 0.81 s.
- `A`  : relative fluorescence increase per single spike (ΔF/F₀).
         Paper OGB mean: 0.052.
- `g`  : saturation parameter, `g = Δ[Ca]_T / ([Ca]_0 + K_d)`. Inverse of
         the spike count at half dye saturation. `g = 0` disables saturation
         (linear convolution limit).
- `F0` : resting-state baseline fluorescence. Arbitrary units; `1.0` gives
         `F` directly as a ΔF/F-ready ratio.
- `η`  : baseline Brownian drift amplitude. `η = 0` gives a flat baseline.
- `σ`  : white measurement-noise standard deviation.
"""
CaModel

@kwdef struct CaModel 
    τ::Float32  = 0.81s
    τr::Float32  = 0.0s
    A::Float32  = 0.052
    g::Float32  = 0.0
    F0::Float32 = 1.0
    η::Float32  = 0.0
    σ::Float32  = 0.0
end

"""
    bin_spikes(spiketimes, dt, interval) -> Vector{Int}

Discretize a spike train onto a regular time grid of step `dt` over
`interval = (t0, t1)`. Returns per-bin spike counts. Spikes outside the
interval are ignored. Implements the Dirac-sum input `s(t) = Σᵢ δ(t - tᵢ)`
integrated over each bin.
"""
function bin_spikes(spiketimes::AbstractVector, dt::Real, interval)
    t0, t1 = interval[1], interval[end]
    n = Int(round((t1 - t0) / dt))
    counts = zeros(Int, n)
    @inbounds for ts in spiketimes
        (ts < t0 || ts >= t1) && continue
        k = Int(floor((ts - t0) / dt)) + 1
        counts[k] += 1
    end
    return counts
end

"""
    calcium_dynamics(spikes_binned, dt, τ, τr) -> Vector{Float64}

Integrate the normalized calcium ODE using exact exponential decay between
bins. When `τr > 0`, calcium rises through an intermediate variable `crise`
with rise time `τr` before decaying with time constant `τ`:

    crise[i] = crise[i-1] * exp(-dt/τr) + spikes[i-1]
    c[i]     = c[i-1]     * exp(-dt/τ)  + crise[i-1]

When `τr = 0` the rise is instantaneous: each spike adds 1 directly to `c`,
recovering the Deneux et al. 2016 single-exponential normalization (one AP →
`c = 1` at rest).
"""
function calcium_dynamics(spikes_binned::T, dt::R, τ::R, τr::R) where {T<:AbstractVector, R<:Real}
    n = length(spikes_binned)
    c = zeros(Float64, n)
    crise = zeros(Float64, n)
    decay = exp(-dt / τ)
    rise = τr > 0 ? exp(-dt / τr) : 0.f0
    α = τr > 0 ? α_doubleexp(τr, τ) : 1.0
    g = τr > 0 ? norm_doubleexp(τr, τ) : 1.0

    @inbounds @simd for i in 2:n
        c[i] = c[i-1] * decay + crise[i-1]
        crise[i] = crise[i-1] * rise + spikes_binned[i-1] * α 
    end
    return g .* c
end

function norm_doubleexp(τr, τd)
    t_p = τr * τd / (τd - τr) * log(τd / τr)
    return 1 / (-exp(-t_p / τr) + exp(-t_p / τd))
end

function α_doubleexp(τr, τd)
    return (τd - τr) / (τd * τr)
end

"""
    baseline_drift(n, dt, η, F0; rng) -> Vector{Float64}

Generate a discretized Brownian-drift baseline `B(t)` of length `n` starting
at `F0`, with increment `η * sqrt(dt) * randn()` per step. When `η = 0` the
baseline is constant at `F0`.
"""
function baseline_drift(n::Int, dt::R, η::R, F0::R; rng=Random.GLOBAL_RNG) where {R<:Real}
    B = fill(float(F0), n)
    η == 0 && return B
    s = η * sqrt(dt)
    random_number = randn(rng, n)
    @inbounds for i in 2:n
        B[i] = B[i-1] + s * random_number[i]
    end
    return B
end

"""
    fluorescence(c, B, A, g, σ; rng) -> Vector{Float64}

Apply the MLspike observation equation
`F = B * (1 + A * c / (1 + g * c)) + σ ε`, combining the saturating
indicator nonlinearity with baseline and additive Gaussian measurement
noise. Reduces to linear scaling when `g = 0`.
"""
function fluorescence(c::AbstractVector, B::AbstractVector, A::Real, g::Real, σ::Real; rng=Random.GLOBAL_RNG)
    n = length(c)
    F = similar(c)
    random_number = randn(rng, n)
    @inbounds for i in 1:n
        nl = A * c[i] / (1 + g * c[i])
        F[i] = B[i] * (1 + nl) + σ * random_number[i]
    end
    return F
end

"""
    calcium_trace(spiketimes, sampling_rate, interval; params=Ca_params, rng)

Simulate a fluorescence trace from a spike train using the Deneux et al.
2016 forward model. `interval = (t0, t1)` in seconds, `sampling_rate` in Hz.
`params` is a NamedTuple of the form [`Ca_params`](@ref).

Pipeline:
1. Bin spikes onto the `dt = 1/sampling_rate` grid.
2. Integrate normalized calcium `c(t)` (exact exponential decay).
3. Generate baseline drift `B(t)` (Brownian, or flat if `η = 0`).
4. Apply the saturating observation equation to obtain `F(t)`.

Returns a NamedTuple `(t, F, c, B)` with the time axis, fluorescence
trace, normalized calcium, and baseline.
"""
function calcium_trace(spiketimes::AbstractVector, sampling_rate::Real, interval::Tuple;
                      params = Ca_params, rng = Random.GLOBAL_RNG)
    dt = 1ms  ## dt for biophysical Calcium dynamics, finer than the output sampling rate to avoid aliasing
    binned = bin_spikes(spiketimes, dt, interval)
    c = calcium_dynamics(binned, dt, params.τ, params.τr)
    B = baseline_drift(length(c), dt, params.η, params.F0; rng=rng)
    F = fluorescence(c, B, params.A, params.g, params.σ; rng=rng)
    t = dt:dt:interval[end]
    c = Itp.scale(Itp.interpolate(c, Itp.BSpline(Itp.Linear())), t)(t[1]:1/sampling_rate:t[end])
    F = Itp.scale(Itp.interpolate(F, Itp.BSpline(Itp.Linear())), t)(t[1]:1/sampling_rate:t[end])
    t = t[1]:1/sampling_rate:t[end]
    return (t=t, F=F, c=c, B=B)
end

export calcium_trace, CaModel, bin_spikes, calcium_dynamics, baseline_drift, fluorescence
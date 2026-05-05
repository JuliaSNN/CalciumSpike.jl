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

Model equations (Deneux et al. 2016, eq. 8, extended with Hill nonlinearity):

    dc/dt = s(t) - c(t) / τ
    dB/dt = η dW(t)
    F(t)  = B(t) * (1 + A * (c0+c)^n / (1 + g*(c0+c)^n) - A*c0^n/(1+g*c0^n)) + σ ε(t)

where `c` is normalized intracellular calcium (c = 0 at rest, c = 1 after a
single AP), `B` is the drifting baseline (Brownian), and `F` the measured
fluorescence. The parameter `c0` is the resting calcium offset (normalized
units) and `n` is the Hill cooperativity coefficient. Setting `n = 1` and
`c0 = 0` recovers the original Michaelis-Menten form. Setting `n > 1` and
`c0 > 0` produces history-dependent fluorescence: early spikes (low baseline
calcium) fall on the shallow left limb of the sigmoid, while accumulated
calcium shifts into the steep central region, amplifying later spikes.
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
- `n`  : Hill cooperativity coefficient for the fluorescence-calcium sigmoid.
         `n = 1` recovers Michaelis-Menten (no cooperativity, current default).
         `n > 1` (typical GCaMP: 2–3) produces a sigmoidal F-Ca curve with
         history-dependent spike reporting: early spikes on the shallow left
         limb produce less ΔF than later spikes in the steep central region.
- `c0` : resting calcium offset in normalized units (same scale as `c`).
         Sets where on the sigmoid curve the indicator sits at rest.
         `c0 = 0` means resting at the foot of the curve (no history effect).
         Increasing `c0` shifts the operating point rightward into the steep
         region. `n = 1, c0 = 0` exactly recovers the original model.
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
    n::Float32  = 1.0
    c0::Float32 = 0.0
    dt::Float32 = 1ms
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
function calcium_dynamics(spikes_binned::T, params::CaModel) where {T<:AbstractVector}
    @unpack τ, τr, dt = params
    n = length(spikes_binned)
    c = zeros(Float32, n)
    crise = zeros(Float32, n)
    decay = exp(-dt / τ) |> Float32
    rise = τr > 0 ? exp(-dt / τr) : 0.f0 |> Float32
    α = τr > 0 ? α_doubleexp(τr, τ) : 1.0f0 |> Float32
    g = τr > 0 ? norm_doubleexp(τr, τ) : 1.0f0 |> Float32
    @inbounds @simd for i in 2:n
        c[i] = c[i-1] * decay + crise[i-1]
        crise[i] = crise[i-1] * rise + spikes_binned[i-1] * α 
    end
    return g .* c
end

function norm_doubleexp(τr, τd)
    t_p = τr * τd / (τd - τr) * log(τd / τr)
    return 1 / (-exp(-t_p / τr) + exp(-t_p / τd)) |> Float32
end

function α_doubleexp(τr, τd)
    return (τd - τr) / (τd * τr) |> Float32
end

"""
    baseline_drift(n, dt, η, F0; rng) -> Vector{Float64}

Generate a discretized Brownian-drift baseline `B(t)` of length `n` starting
at `F0`, with increment `η * sqrt(dt) * randn()` per step. When `η = 0` the
baseline is constant at `F0`.
"""
function baseline_drift(n::Int, params::CaModel; rng=Random.GLOBAL_RNG)
    @unpack dt, η, F0 = params
    B = fill(float(F0), n)
    η == 0 && return B
    s = η * sqrt(dt)
    random_number = randn(rng, n)
    @turbo for i in 2:n
        B[i] = B[i-1] + s * random_number[i]
    end
    return B
end

"""
    fluorescence(c, B, params; rng) -> Vector{Float32}

Apply the observation equation combining the Hill-sigmoid indicator
nonlinearity with baseline and additive Gaussian measurement noise:

    F = B * (1 + A*(c0+c)^n/(1+g*(c0+c)^n) - A*c0^n/(1+g*c0^n)) + σ ε

The subtracted rest term normalizes so that F(c=0) = B·(1+0) = B.
Setting `n=1, c0=0` recovers the original Michaelis-Menten equation.
Setting `n>1` introduces sigmoidal cooperativity; combining with `c0>0`
places the resting operating point on the rising limb so that accumulated
calcium amplifies the fluorescence response to subsequent spikes
(history dependence, as in GCaMP6f; see Demas et al. 2021).
"""
function fluorescence(c::T, B::T, params::CaModel; rng=Random.GLOBAL_RNG) where {T<:AbstractVector}
    @unpack A, g, σ, n, c0 = params
    nn = length(c)
    F = similar(c)
    random_number = randn(rng, nn)
    F_rest = A * c0^n / (1 + g * c0^n)
    @inbounds for i in 1:nn
        c_eff = c0 + c[i]
        nl = A * c_eff^n / (1 + g * c_eff^n) - F_rest
        F[i] = B[i] * (1 + nl) + σ * random_number[i]
    end
    return F
end

"""
    calcium_trace(spiketimes::Vector{Float32}, sampling_rate::Real, interval; params=Ca_params, rng=Random.GLOBAL_RNG)

Simulate a fluorescence trace from a spike train using the Deneux et al.
2016 forward model. `interval = (t0, t1)` in seconds, `sampling_rate` in Hz.
`params` is a [`CaModel`](@ref) instance.

Pipeline:
1. Bin spikes onto the `dt = 1/sampling_rate` grid.
2. Integrate normalized calcium `c(t)` (exact exponential decay).
3. Generate baseline drift `B(t)` (Brownian, or flat if `η = 0`).
4. Apply the saturating observation equation to obtain `F(t)`.

Returns a NamedTuple `(t, F, c, B)` with the time axis, fluorescence
trace, normalized calcium, and baseline.
"""
function calcium_trace(spiketimes::Vector{Float32}, sampling_rate::Real, interval;
                      params = Ca_params, rng = Random.GLOBAL_RNG)
    dt = 1ms  ## dt for biophysical Calcium dynamics, finer than the output sampling rate to avoid aliasing
    binned = bin_spikes(spiketimes, dt, interval)
    c = calcium_dynamics(binned, params)
    B = baseline_drift(length(c), params; rng=rng)
    F = fluorescence(c, B, params; rng=rng)
    t = dt:dt:interval[end]
    c = Itp.scale(Itp.interpolate(c, Itp.BSpline(Itp.Linear())), t)(t[1]:1/sampling_rate:t[end])
    F = Itp.scale(Itp.interpolate(F, Itp.BSpline(Itp.Linear())), t)(t[1]:1/sampling_rate:t[end])
    t = t[1]:1/sampling_rate:t[end]
    return (t=t, F=F, c=c, B=B)
end

"""
    calcium_trace(spiketimes::Vector{Vector{Float32}}, sampling_rate::Real, interval; params=Ca_params, rng=Random.GLOBAL_RNG)

    Simulate a fluorescence trace from a vector of spike trains using the Deneux et al.
    2016 forward model. `interval = (t0, t1)` in seconds, `sampling_rate` in Hz.
`params` is a [`CaModel`](@ref) instance.

Pipeline:
1. Bin spikes onto the `dt = 1/sampling_rate` grid.
2. Integrate normalized calcium `c(t)` (exact exponential decay).
3. Generate baseline drift `B(t)` (Brownian, or flat if `η = 0`).
4. Apply the saturating observation equation to obtain `F(t)`.

Returns a tuple of vectors `(Fs, t)` with the fluorescence traces and time axis.
"""
function calcium_trace(spiketimes::Vector{Vector{Float32}}, sampling_rate::Real, interval;
        params = Ca_params, rng = Random.GLOBAL_RNG)
    dt = 1ms  ## dt for biophysical Calcium dynamics, finer than the output sampling rate to avoid aliasing
    t = 0:dt:interval[end]
    tsr = t[1]:1/sampling_rate:t[end]
    Fs = map(spiketimes) do spiketime
        binned = bin_spikes(spiketime, dt, interval)
        c = calcium_dynamics(binned, params)
        B = baseline_drift(length(c), params; rng=rng)
        F = fluorescence(c, B, params; rng=rng)
        t0 = range(t[1], stop=t[end], length=length(F)) 
        F = Itp.scale(Itp.interpolate(F, Itp.BSpline(Itp.Linear())), t0)(tsr)
    end
    return Fs, tsr
end

"""
    gcamp6_kernel(;τr, τd, interval) -> Vector{Float32}

Generate a discretized gCaMP6s kernel of the form `g(t) = (exp(-t/τd) - exp(-t/τr))`
normalized to unit integral. `interval` is used to determine the kernel time axis and bin width. The kernel is truncated at `kernel_length = (τr + τd) * 10` to ensure it captures the full decay while avoiding unnecessary computation of negligible tails.
"""
function gcamp6_kernel(;τr=100ms, τd=2s, interval, kwargs...)
    kernel_length = (τr+τd) * 10  # Length of the kernel in ms
    bin_width = step(interval) # kernel window in ms
    kernel_time = Float32.(0.0:bin_width:kernel_length)
    τr = Float32(τr)
    τd = Float32(τd)
    gcamp6_kernel = [(exp(-t / τd)) - exp(-t / τr) for t in kernel_time]
    gcamp6_kernel ./= bin_width * sum(gcamp6_kernel)
    return gcamp6_kernel
end



export calcium_trace, CaModel, bin_spikes, calcium_dynamics, baseline_drift, fluorescence, gcamp6_kernel
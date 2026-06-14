"""
    Ca_model

Forward model of calcium fluorescence from spike trains, following the
exact physiological model of Deneux et al. 2016 (equation 17), as used by
MLspike.

Reference:
    Deneux, T., Kaszas, A., Szalay, G., Katona, G., Lakner, T., Grinvald, A.,
    Rózsa, B., & Vanzetta, I. (2016). Accurate spike estimation from noisy
    calcium signals for ultrafast three-dimensional imaging of large
    neuronal populations in vivo. *Nature Communications*, 7, 12190.
    doi:10.1038/ncomms12190

Model equations (Deneux et al. 2016, eq. 17):

    dc/dt = -c(t) / τ + s(t)
    dp/dt = (1/τr) * ( ((c0+c)^n - c0^n) / (1 + g*((c0+c)^n - c0^n)) - p )
    dB/dt = ζ dW(t)
    F(t)  = B(t) * (1 + A * p(t)) + σ ε(t)

where `c` is normalized intracellular calcium (c = 0 at rest, c = 1 after a
single AP), `p` is the fraction of indicator bound to calcium (a filtered,
Hill-transformed version of the calcium increment), `B` is the drifting
baseline (Brownian), and `F` the measured fluorescence. The parameter `c0`
is the resting calcium offset (normalized units) and `n` is the Hill
cooperativity coefficient. Setting `n = 1` and `c0 = 0` and `τr = 0`
recovers the original Michaelis-Menten linear limit
`F = B * (1 + A * c / (1 + g*c)) + σε`. Setting `n > 1` and `c0 > 0`
produces history-dependent fluorescence: early spikes (low baseline
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
- `τr`: indicator binding time constant (Deneux `t_on`, seconds). Sets
         the rise time of the bound-indicator fraction `p` toward its
         instantaneous Hill target `f(c)`. `τr = 0` collapses `p` to
         `f(c)` evaluated pointwise (instantaneous binding, no ODE);
         this is the default and recovers the original Deneux limit when
         combined with `n = 1, c0 = 0`.
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
         region. `n = 1, c0 = 0, τr = 0` exactly recovers the original model.
"""
CaModel

@kwdef struct CaModel
    τ::Float32  = 0.81s
    τr::Float32 = 0.0s
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
    calcium_dynamics(spikes_binned, params) -> Vector{Float32}

Integrate the normalized calcium ODE `dc/dt = -c/τ + s(t)` using exact
exponential decay between bins:

    c[i] = c[i-1] * exp(-dt/τ) + spikes[i-1]

This is the pure single-exponential form of the Deneux et al. 2016 model
(eq. 17, calcium line). The indicator binding rise time `τr` is handled
by the separate state variable `p` in [`indicator_dynamics`](@ref); it
does NOT appear here.
"""
function calcium_dynamics(spikes_binned::T, params::CaModel) where {T<:AbstractVector}
    @unpack τ, dt = params
    n = length(spikes_binned)
    c = zeros(Float32, n)
    decay = exp(-dt / τ) |> Float32
    @inbounds @simd for i in 2:n
        c[i] = c[i-1] * decay + spikes_binned[i-1]
    end
    return c
end

"""
    indicator_dynamics(c, params) -> Vector{Float32}

Compute the fraction of indicator bound to calcium `p(t)` from the
normalized calcium trace `c(t)`, following Deneux et al. 2016 eq. 17:

    dp/dt = (1/τr) * ( f(c) - p )

with the Hill-transformed calcium *increment* (note: increment, not raw
power)

    f(c) = ((c0 + c)^n - c0^n) / (1 + g * ((c0 + c)^n - c0^n))

The ODE is integrated with the exact-step formula for piecewise-constant
input over each bin of width `dt`:

    p[i] = p[i-1] * exp(-dt/τr) + (1 - exp(-dt/τr)) * f(c[i-1])

When `τr ≤ 0` the binding is treated as instantaneous and the ODE
collapses to the algebraic limit `p[i] = f(c[i])` evaluated at each
sample (no filtering).

Limit check: with `n = 1, c0 = 0, τr = 0`, we have
`f(c) = c / (1 + g*c)` (since `(0+c)^1 - 0^1 = c`), and `p = f(c)`
pointwise, so the observation equation in [`fluorescence`](@ref) reduces
to `F = B * (1 + A * c / (1 + g*c))`, the original Michaelis-Menten form.
"""
function indicator_dynamics(c::T, params::CaModel) where {T<:AbstractVector}
    @unpack A, g, n, c0, τr, dt = params
    N = length(c)
    p = zeros(Float32, N)

    # Hill-transformed calcium increment (Deneux numerator)
    function f(ci)
        x = (c0 + ci)^n - c0^n
        x / (1 + g * x)
    end

    if τr <= 0
        # instantaneous binding: p = f(c) at each step
        @inbounds for i in 1:N
            p[i] = f(c[i])
        end
    else
        decay = exp(-dt / τr)
        @inbounds for i in 2:N
            p[i] = p[i-1] * decay + (1 - decay) * f(c[i-1])
        end
    end
    return p
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

Apply the Deneux et al. 2016 (eq. 17) observation equation. The Hill
nonlinearity is *not* applied here; it is already absorbed into the
indicator-binding state `p(t)`, which is computed internally by
[`indicator_dynamics`](@ref). The observation is therefore linear in `p`:

    F = B * (1 + A * p) + σ ε

where `p = indicator_dynamics(c, params)`.

Limit check: with `τr = 0, n = 1, c0 = 0`, `indicator_dynamics`
returns `p = c / (1 + g*c)` pointwise, so

    F = B * (1 + A * c / (1 + g*c)) + σ ε,

the original Michaelis-Menten form. Setting `n > 1` and `c0 > 0`
introduces sigmoidal cooperativity and history dependence (the operating
point sits on the rising limb so accumulated calcium amplifies subsequent
ΔF; see Demas et al. 2021). Setting `τr > 0` adds an indicator-binding
rise time, low-pass filtering `f(c)` with time constant `τr`.
"""
function fluorescence(c::T, B::T, params::CaModel; rng=Random.GLOBAL_RNG) where {T<:AbstractVector}
    @unpack A, σ = params
    nn = length(c)
    p = indicator_dynamics(c, params)
    F = similar(c)
    random_number = randn(rng, nn)
    @inbounds for i in 1:nn
        F[i] = B[i] * (1 + A * p[i]) + σ * random_number[i]
    end
    return F
end

"""
    calcium_trace(spiketimes::Vector{Float32}, sampling_rate::Real, interval; params=Ca_params, rng=Random.GLOBAL_RNG)

Simulate a fluorescence trace from a spike train using the Deneux et al.
2016 forward model (eq. 17). `interval = (t0, t1)` in seconds,
`sampling_rate` in Hz. `params` is a [`CaModel`](@ref) instance.

Pipeline:
1. Bin spikes onto the `dt = 1/sampling_rate` grid.
2. Integrate normalized calcium `c(t)` (exact single-exponential decay).
3. Generate baseline drift `B(t)` (Brownian, or flat if `η = 0`).
4. Apply the Deneux observation equation `F = B * (1 + A * p) + σ ε`,
   where `p` (indicator-bound fraction) is computed internally by
   [`indicator_dynamics`](@ref) from `c`.

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
    2016 forward model (eq. 17). `interval = (t0, t1)` in seconds, `sampling_rate` in Hz.
`params` is a [`CaModel`](@ref) instance.

Pipeline:
1. Bin spikes onto the `dt = 1/sampling_rate` grid.
2. Integrate normalized calcium `c(t)` (exact single-exponential decay).
3. Generate baseline drift `B(t)` (Brownian, or flat if `η = 0`).
4. Apply the Deneux observation equation `F = B * (1 + A * p) + σ ε`,
   where `p` (indicator-bound fraction) is computed internally by
   [`indicator_dynamics`](@ref) from `c`.

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



export calcium_trace, CaModel, bin_spikes, calcium_dynamics, indicator_dynamics, baseline_drift, fluorescence, gcamp6_kernel

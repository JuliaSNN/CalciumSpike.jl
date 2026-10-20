# Biophysical Forward Model

Implementation of the calcium fluorescence forward model from
Deneux et al. 2016 (Nat Commun 7:12190).

```@meta
CurrentModule = CalciumSpike
```

## Model equations

Normalized intracellular calcium ``c(t)`` is driven by the spike input
``s(t) = \sum_i \delta(t - t_i)`` and decays with time constant ``\tau``:

```math
\frac{dc}{dt} = s(t) - \frac{c(t)}{\tau}
```

When a non-zero rise time ``\tau_r`` is specified, calcium follows a
dual-exponential kernel through an intermediate rise variable ``c_\text{rise}``:

```math
\frac{dc_\text{rise}}{dt} = s(t) - \frac{c_\text{rise}}{\tau_r}, \qquad
\frac{dc}{dt} = c_\text{rise} - \frac{c}{\tau}
```

The baseline ``B(t)`` is a Brownian drift:

```math
dB = \eta \, dW(t)
```

Measured fluorescence combines the saturating indicator nonlinearity with
additive Gaussian noise:

```math
F(t) = B(t)\!\left(1 + \frac{A \, c(t)}{1 + g \, c(t)}\right) + \sigma\,\varepsilon(t)
```

Setting ``g = 0`` and ``\eta = 0`` recovers the linear convolution limit.

## Randomness and reproducibility of the forward model

By default `calcium_trace` is random: the baseline drift (`η > 0`) and the
measurement noise (`σ > 0`) come from the global RNG, so each call gives a new
noise realisation (`Random.seed!` before the call makes it repeatable).

To opt in to bit-for-bit reproducibility pass an integer `seed`:

```julia
Fs, t = calcium_trace(spikes, 50.0, (0.0, 30.0); params, seed = 42)
```

The call builds a local `Random.Xoshiro(42)` and uses it for all of its noise,
so the same spikes, `params` and `seed` give identical output regardless of the
global RNG state, earlier calls, or the number of threads. For a
`Vector{Vector{Float32}}` input the neurons are processed serially from that one
stream; within one seed every call with the same spike input gives the same
noise realisation (common random numbers across perturbation conditions).
`seed` and `rng` are mutually exclusive. `biophysical_calcium(...; seed)` seeds
the calcium noise only, not the Poisson spike trains, and
`generate_synthetic_data(...; seed)` seeds the synthetic data.

## API

```@autodocs
Modules = [CalciumSpike]
Order   = [:type, :function]
Filter  = t -> nameof(t) in (:CaModel, :bin_spikes, :calcium_dynamics,
                              :baseline_drift, :fluorescence, :calcium_trace)
```

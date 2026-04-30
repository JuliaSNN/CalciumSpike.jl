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

## API

```@autodocs
Modules = [CalciumSpike]
Order   = [:type, :function]
Filter  = t -> nameof(t) in (:CaModel, :bin_spikes, :calcium_dynamics,
                              :baseline_drift, :fluorescence, :calcium_trace)
```

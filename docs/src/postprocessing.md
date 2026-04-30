# Signal Post-Processing

Standard pipeline for converting raw fluorescence to estimated firing rate.

```@meta
CurrentModule = CalciumSpike
```

## Method

ΔF/F is computed from the raw trace ``F(t)`` using a baseline ``F_0``
estimated from the lower quantile of the recording:

```math
\frac{\Delta F}{F}(t) = \frac{F(t) - F_0}{F_0}
```

Analytical deconvolution under the single-exponential assumption recovers an
estimated firing rate ``r(t)`` directly from the ΔF/F signal:

```math
r(t) \approx \frac{1}{A}\!\left(\frac{d(\Delta F/F)}{dt} + \frac{\Delta F/F}{\tau}\right)
```

A Gaussian kernel with width ``\sigma_\text{smooth}`` is applied afterwards to
suppress differentiation noise.

## API

```@autodocs
Modules = [CalciumSpike]
Order   = [:type, :function]
Filter  = t -> nameof(t) in (:CaPostProcess, :delta_f_over_f, :gaussian_smooth,
                              :deconvolve_df_f, :calcium_postprocess)
```

# Signal Post-Processing

Standard pipeline for converting raw fluorescence to estimated firing rate.

```@meta
CurrentModule = CalciumSpike
```

## Method

### Baseline estimation (`delta_f_over_f`)

A causal running median over a window ``W`` (default 3 s) is applied to the raw
trace; the 20th percentile of the resulting trajectory is taken as a single
scalar baseline ``F_0``:

```math
\tilde{F}(t) = \operatorname{median}_{s \in [t-W,\, t]} F(s)
\qquad
F_0 = Q_{0.20}\!\bigl(\tilde{F}\bigr)
```

ΔF/F follows:

```math
\frac{\Delta F}{F}(t) = \frac{F(t) - F_0}{F_0}
```

This matches the Sophie/Deneux MATLAB pipeline
(`fn_filt` running median → 20th-percentile F0) and is more robust than a
simple mean over a fixed warmup window.

### Deconvolution

Analytical deconvolution under the single-exponential assumption recovers an
estimated firing rate ``r(t)`` directly from the ΔF/F signal:

```math
r(t) \approx \frac{1}{A}\!\left(\frac{d(\Delta F/F)}{dt} + \frac{\Delta F/F}{\tau}\right)
```

### Smoothing (`gaussian_smooth`)

A Gaussian kernel with width ``\sigma_\text{smooth}`` suppresses differentiation noise.
The kernel can be **symmetric** (default), **causal** (`skewed=:left`, uses only past
samples — no lookahead, no leakage), or **anti-causal** (`skewed=:right`).
The causal kernel matches Sophie's `gaussianfilter(..., causal=1)`.

## API

```@autodocs
Modules = [CalciumSpike]
Order   = [:type, :function]
Filter  = t -> nameof(t) in (:CaPostProcess, :delta_f_over_f, :gaussian_smooth,
                              :deconvolve_df_f, :calcium_postprocess)
```

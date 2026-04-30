# Noise Correction

Trial-averaged population correlation analysis with noise-correlation correction.

```@meta
CurrentModule = CalciumSpike
```

## Method

For each neuron pair ``(i, j)``, the raw cross-correlation ``\rho_0`` is
computed across split-half trial groups. Noise correlations
``\rho^\text{noise}_{i}`` and ``\rho^\text{noise}_{j}`` are estimated from
within-condition variability. The corrected correlation is:

```math
\rho_{ij} = \frac{\rho_0}{\sqrt{|\rho^\text{noise}_i \cdot \rho^\text{noise}_j|}}
```

Pairs where the correction factor is below a threshold (< 0.05) are set to
zero to avoid amplifying noise-dominated estimates.

## Functions

```@autodocs
Modules = [CalciumSpike]
Order   = [:function]
Filter  = t -> t in (split_mean, generate_synthetic_data, noise_correction)
```

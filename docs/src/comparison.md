# Model Comparison

Population-level simulation routines for benchmarking spike-inference approaches.

```@meta
CurrentModule = CalciumSpike
```

## Workflow

```
spike train
    ├─── CaModel (Deneux 2016) ──► ΔF/F ──► deconvolution ──► Pearson r
    └─── double-exp gCaMP6    ──► ΔF/F ──► deconvolution ──► Pearson r
                                                                    │
                                                         ground-truth firing rate
```

Both approaches are compared against the ground-truth firing rate using Pearson
correlation. The `β` parameter of the inhomogeneous Poisson model controls
firing irregularity; use [`isi_beta`](@ref) to measure the resulting ISI-CV.

## Functions

```@autodocs
Modules = [CalciumSpike]
Order   = [:function]
Filter  = t -> t in (ca_fr_correlation, run_comparison, biophysical_calcium, evaluate_deconvolution, isi_beta)
```

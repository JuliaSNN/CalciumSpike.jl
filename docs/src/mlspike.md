# MLSpike Integration

Julia interface to the MATLAB MLSpike toolbox (Deneux et al. 2016).

```@meta
CurrentModule = CalciumSpike
```

!!! warning "Lazy loading"
    MLSpike functions are loaded only when `MATLAB.jl` is present in the
    environment. Load MATLAB.jl before CalciumSpike, or trigger loading with
    `import MATLAB` after `using CalciumSpike`.

## Requirements

- MATLAB version ≤ 2022
- Optimization Toolbox
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- [spikes](https://github.com/MLspike/spikes) library
- [brick](https://github.com/thomasdeneux/brick) library

## Setup

```julia
import MATLAB
using CalciumSpike

activate_MLSpike("/path/to/mlspike_root/")
```

The root directory must contain `spikes/` and `brick/` subdirectories.

## Pipeline

```
autocalibration(calcium)          # estimate a, tau, sigma from data
    └─► estimation_params(...)    # build MATLAB parameter struct
            └─► spike_estimate(calcium, par)   # infer spike times
```

Or use the high-level wrapper:

```julia
result = MLspike_estimate(calcium)
# result.spikest, result.fit, result.drift
```

## Functions

```@autodocs
Modules = [CalciumSpike]
Order   = [:function]
Filter  = t -> nameof(t) in (:activate_MLSpike, :generate_spike_train, :spk_calcium,
                              :autocalibration, :estimation_params, :spike_estimate,
                              :MLspike_estimate, :evaluate_MLspike)
```

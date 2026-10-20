# CalciumSpike changelog

## 1.1.0

New opt-in feature; no change for code that does not use it.

- `calcium_trace` (both the `Vector{Float32}` and `Vector{Vector{Float32}}` methods) takes `seed = nothing`. The default
  stays random: noise comes from the global RNG, fresh at every call, as in 1.0.0. With `seed::Integer` the call builds a
  local `Random.Xoshiro(seed)` and draws all its noise from it, so equal spikes, `params` and `seed` give bitwise-identical
  output independent of the global RNG state, earlier calls and thread count. The vector method draws neuron by neuron,
  serially, from that single stream (common random numbers across conditions that share the spike input). `rng` and `seed`
  together is an `ArgumentError`. The `rng` default is now `nothing` (meaning the global RNG) instead of `Random.GLOBAL_RNG`.
- `generate_synthetic_data` takes `rng` and `seed` (default: global RNG, as before).
- `biophysical_calcium` takes `seed`: seeds the calcium noise per neuron (derived from a local stream, thread-independent).
  The Poisson spike trains of the `SNNModels` simulation are not seeded by it.
- First test suite (`test/runtests.jl`) covering the seed behaviour.
- `CaModel` is unchanged: the seed is a keyword of the calls, not a model parameter, so positional constructors and saved
  parameters keep working.

## 1.0.0

Initial release.

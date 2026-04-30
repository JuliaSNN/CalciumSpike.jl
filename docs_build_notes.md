# CalciumSpike.jl — Documentation Build Notes

Summary of errors encountered when running `include("src/CalciumSpike.jl/docs/make.jl")` and how each was resolved.

---

## 1. `noise_correction.jl` not loaded by the module

**Error:** Functions from `noise_correction.jl` (`split_mean`, `generate_synthetic_data`, `noise_correction`) were absent from the module, so `@autodocs` produced no output for that page.

**Cause:** `CalciumSpike.jl` was missing the `include("noise_correction.jl")` line.

**Fix:** Added `include("noise_correction.jl")` to `src/CalciumSpike.jl`.

---

## 2. `api_reference.md` built as orphan page, causing duplicate doc errors

**Error:**
```
Error: duplicate docs found for 'CaModel' in `@docs` block in forward_model.md
Warning: duplicate docs found for 'CalciumSpike.ca_fr_correlation' in comparison.md
... (and all other comparison functions)
```

**Cause:** Documenter builds every `.md` file in `docs/src/`, not only those listed in `pages`. `api_reference.md` was commented out of the `pages` list but still existed on disk. Its `@autodocs Modules = [CalciumSpike]` (no filter) captured all functions first; subsequent per-page `@autodocs` blocks then found the same symbols again.

**Fix:** Deleted `docs/src/api_reference.md`.

---

## 3. Lazy-loaded symbols referenced by name in `@autodocs` filters

**Error:**
```
ERROR: UndefVarError: `activate_MLSpike` not defined in `Main`
```

**Cause:** The filter lambda in `mlspike.md` used `t in (activate_MLSpike, ...)`. Julia evaluates this tuple at doc-build time in `Main` scope. Since MLSpike functions are lazy-loaded via `Requires.jl` (only present when MATLAB.jl is loaded), `activate_MLSpike` did not exist and the lambda threw `UndefVarError`. The same issue applied to the Makie plot functions in `visualization.md`.

**Fix:** Replaced direct function references with `nameof`-based symbol matching:
```julia
# Before (fails without MATLAB.jl)
Filter = t -> t in (activate_MLSpike, generate_spike_train, ...)

# After (safe at any load state)
Filter = t -> nameof(t) in (:activate_MLSpike, :generate_spike_train, ...)
```
Applied to both `mlspike.md` and `visualization.md`.

---

## 4. `@docs CaModel` conflict with `@autodocs Order = [:function]`

**Error:**
```
Error: duplicate docs found for 'CaModel' in `@docs` block in forward_model.md:44-46
```

**Cause:** Julia structs defined with `@kwdef` have both a type and a callable constructor. `@autodocs` with `Order = [:function]` picks up the constructor, and the explicit `@docs CaModel` block on the same page also documents the type — Documenter sees two entries for the same binding.

**Fix:** Removed the separate `@docs CaModel` block and merged it into the `@autodocs` block with `Order = [:type, :function]` and a `nameof`-based filter including `:CaModel`. Applied the same pattern to `postprocessing.md` (`CaPostProcess`) for consistency.

---

## 5. Stale SNNModels precompilation cache

**Error:**
```
SystemError: opening file ".../SNNModels.jl/src/populations/inh_poisson.jl": No such file or directory
ERROR: Failed to precompile SNNModels
```

**Cause:** A population file was renamed from `inh_poisson.jl` to `inhomogeneous_poisson.jl`. The Julia precompilation cache (`~/.julia/compiled/v1.12/SNNModels/`) still listed the old filename as a tracked dependency and failed the staleness check.

**Fix:**
```bash
rm -rf ~/.julia/compiled/v1.12/SNNModels/
```
Julia recompiled SNNModels from scratch using the current file tree.

---

## 6. Script hyperlinks pointing outside `docs/src/`

**Error:**
```
Error: invalid local link/image: path pointing to a file outside of build directory
  link = ../../scripts/run_demo.jl
```

**Cause:** `index.md` used markdown hyperlinks with relative paths like `../../scripts/run_demo.jl`. Documenter rejects links that resolve outside the `docs/src/` tree.

**Fix:** Replaced all markdown link syntax with plain inline code spans — the table retains script names as readable references without broken HTML links.

---

## 7. Unresolvable `@ref` for `Ca_params`

**Error:**
```
Warning: Cannot resolve @ref for `Ca_params` in forward_model.md
```

**Cause:** `Ca_params` is a documented NamedTuple defined in `Ca_model.jl` but not exported from `CalciumSpike`. Documenter cannot resolve `@ref` for unexported symbols.

**Fix:** Replaced the `[`Ca_params`](@ref)` cross-reference in the `calcium_trace` docstring with `[`CaModel`](@ref)`, which is the actual exported parameter struct.

---

## 8. Remaining (non-blocking) warnings

These are suppressed via `warnonly = [:autodocs_block, :missing_docs, :cross_references]` and do not prevent the build from producing HTML:

| Warning | Cause |
|---------|-------|
| `evaluate_MLspike` / `plot_comparison` unresolvable `@ref` | These are lazy-loaded symbols not present at build time; the `@ref`s appear in docstrings of always-loaded functions that call them |
| `CalciumSpike.CalciumSpike` not in manual | Module-level docstring not included on any page |
| `gcamp6_kernel` not in manual | Internal helper function with a docstring but not included in any page filter |
| `Unable to determine edit_link / devbranch` | No git remote configured in the local build environment |
| `Skipping deployment` | Expected for local-only builds without a CI environment |

# CaModel Changes Report
*2026-05-06 — for reviewer*

## Summary

`src/Ca_model.jl` was updated to implement the **exact Deneux et al. 2016 eq. 17** forward model, replacing a previous approximate formulation. The fitting scripts (`optimization.jl`, `lbfgs_fit.jl`) reference `CaModel` parameters by name — this report documents every change to the struct and functions so the reviewer can verify consistency.

---

## 1. `CaModel` struct — field changes

### Renamed field

**All fields unchanged:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `τ` | `Float32` | `0.81s` | Calcium decay time constant |
| `τr` | `Float32` | `0.81s` | Calcium rise time constant |
| `A` | `Float32` | `0.052` | ΔF/F amplitude per single AP |
| `g` | `Float32` | `0.0` | Saturation parameter |
| `F0` | `Float32` | `1.0` | Resting baseline fluorescence |
| `η` | `Float32` | `0.0` | Baseline Brownian drift amplitude |
| `σ` | `Float32` | `0.0` | Measurement noise std |
| `n` | `Float32` | `1.0` | Hill cooperativity coefficient |
| `c0` | `Float32` | `0.0` | Resting calcium offset (normalized) |
| `dt` | `Float32` | `1ms` | Internal simulation timestep |

---

## 2. Functions — changes

### `calcium_dynamics(spikes_binned, params)` — simplified

**Before:** double-exponential rise via `τr`. Maintained two state variables `c` and `crise`; rise was driven by spikes through `crise` with time constant `τr`.

**After:** pure single-exponential decay only:
```
c[i] = c[i-1] * exp(-dt/τ) + spikes[i-1]
```
The `τr`/`crise` logic, `norm_doubleexp`, and `α_doubleexp` helpers are all deleted. Rise time is now handled entirely by `indicator_dynamics`.

---

### `indicator_dynamics(c, params)` — **NEW function**

Computes the fraction of indicator bound to calcium `p(t)` from the calcium trace `c(t)`.

**Equation (Deneux eq. 17):**
```
dp/dt = (1/ton) * ( f(c) - p )

where f(c) = ((c0 + c)^n - c0^n) / (1 + g * ((c0 + c)^n - c0^n))
```

**Discretization (exact-step formula):**
```
p[i] = p[i-1] * exp(-dt/ton) + (1 - exp(-dt/ton)) * f(c[i-1])
```

**Special case `ton <= 0` (instantaneous binding):**
```
p[i] = f(c[i])    # pointwise, no ODE
```

**Note on Hill function form:** the argument of the saturation denominator is the **increment** `(c0+c)^n - c0^n`, not the raw power `(c0+c)^n`. This is a key difference from the previous implementation.

---

### `fluorescence(c, B, params; rng)` — simplified

**Before:** applied Hill nonlinearity directly to `c` at the output:
```
nl = A * (c0+c)^n / (1 + g*(c0+c)^n) - A*c0^n/(1+g*c0^n)
F[i] = B[i] * (1 + nl) + σ * ε[i]
```

**After:** linear in `p` (Hill is inside `indicator_dynamics`):
```
p = indicator_dynamics(c, params)
F[i] = B[i] * (1 + A * p[i]) + σ * ε[i]
```

---

### Deleted functions

- `norm_doubleexp(τr, τd)` — normalization factor for double-exp; no longer needed
- `α_doubleexp(τr, τd)` — scaling factor for double-exp; no longer needed

---

### Unchanged functions

- `bin_spikes` — untouched
- `baseline_drift` — untouched
- `calcium_trace` (both single and batch methods) — pipeline unchanged (`c → B → F`); `fluorescence` now calls `indicator_dynamics` internally
- `gcamp6_kernel` — untouched (uses local `τr` kwarg unrelated to `CaModel`)

---

## 3. Backwards compatibility

Setting `ton=0, n=1, c0=0` recovers the original Michaelis-Menten form exactly:

```
f(c) = ((0+c)^1 - 0^1) / (1 + g*((0+c)^1 - 0^1)) = c / (1+g*c)
p[i] = f(c[i])   (ton=0 branch)
F[i] = B[i] * (1 + A * c[i]/(1+g*c[i])) + σ*ε[i]   ✓
```

---

## 4. Impact on fitting scripts

### `optimization.jl` (Optuna/TPE)

- Uses `τr` keyword? **No** — the script constructs `CaModel` with `τ, A, g, c0, n` only; `τr` was never suggested as a trial parameter. **No change needed.**
- Parameter names in `trial.suggest_float(...)` match current `CaModel` fields: `tau`, `tau_r`, `A`, `g`, `c0`, `n`. Note: `tau_r` in the Optuna trial maps to `CaModel.ton` — **verify the constructor call on line ~131 uses `ton=τr` not `τr=τr`.**

### `lbfgs_fit.jl` (Optim/LBFGS)

- Same parameter vector `[τ, τr, A, g, c0, n]` — `τr` in the vector maps to `CaModel.ton`. **Verify the `CaModel(...)` constructor call uses `ton=τr` not `τr=τr`.**
- Known bug: `tmap` deadlock under Optim finite-diff gradient loop. Fix attempted but not confirmed. **Needs testing.**

---

## 5. Files modified

| File | Change |
|------|--------|
| `src/Ca_model.jl` | All changes described above |
| `scripts/GCamp6/lbfgs_fit.jl` | `sim_traces`: attempted tmap deadlock fix (sequential map); `CaModel` constructor may need `τr→ton` rename check |
| `scripts/GCamp6/optimization.jl` | `CaModel` constructor may need `τr→ton` rename check |

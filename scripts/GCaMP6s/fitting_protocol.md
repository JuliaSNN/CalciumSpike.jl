# GCaMP6s CaModel Fitting Protocol

## Target data (Chen et al. 2013)

Two complementary constraints:

**1. Half-decay time curve** — time from peak to half-peak as a function of spike count (8 points, `DATA_X`/`T_HALF_EMPIRICAL` in `gcamp6s_data.jl`):

| APs | 1 | 2 | 3 | 5 | 10 | 20 | 40 | 80 |
|-----|---|---|---|---|----|----|----|----|
| t½ (s) | 1.03 | 1.43 | 1.49 | 1.58 | 1.81 | 2.27 | 2.73 | 3.28 |

**2. Transient waveforms** — digitized ΔF/F time courses for a 1 AP and 10 AP burst (21 time points each, −0.4 s to +3.0 s relative to spike onset).

---

## Free parameters

| Parameter | Meaning | Range |
|-----------|---------|-------|
| `τ` | Calcium decay time constant | 0.5 – 3.0 s |
| `τr` | Rise time (double-exponential) | 0.01 – 0.5 s |
| `A` | ΔF/F amplitude per single AP | 0.05 – 1.0 |
| `g` | Saturation (inverse half-sat spike count) | 0.0 – 0.3 |
| `c0` | Resting calcium offset (normalized) | 0.0 – 0.5 |
| `n` | Hill cooperativity coefficient | 1.0 – 4.0 |

Fixed: `σ = 0` (no observation noise), `η = 0` (no baseline drift), `F0 = 1`.

---

## Forward model

For a given parameter set, simulate at 1 ms internal resolution, downsampled to 50 Hz:

$$\frac{dc}{dt} = s(t) - \frac{c}{\tau} \qquad \text{(double-exp rise via } \tau_r\text{)}$$

$$F(t) = F_0 \left(1 + \frac{A\,(c_0+c)^n}{1 + g\,(c_0+c)^n} - \frac{A\,c_0^n}{1+g\,c_0^n}\right)$$

$$\frac{\Delta F}{F} = \frac{F - F_0^\text{base}}{F_0^\text{base}}, \qquad F_0^\text{base} = \langle F(t > 10\,\text{s}) \rangle$$

Spikes are placed as an instantaneous burst at $t = 4$ s; the 30 s post-spike window ensures the transient decays fully before baseline estimation.

---

## Objective (minimize)

Three MSE terms summed with equal weight:

$$\mathcal{L}_\text{half} = \frac{1}{8}\sum_{i=1}^{8}\left(\log \hat{t}_{1/2}^{(i)} - \log t_{1/2}^{(i)}\right)^2$$

$$\mathcal{L}_\text{1ap} = \frac{1}{J}\sum_{j=1}^{J}\left(\hat{F}_\text{norm}(t_j) - F^\text{emp}_\text{norm}(t_j)\right)^2$$

$$\mathcal{L}_\text{10ap} = \text{same as } \mathcal{L}_\text{1ap} \text{ for the 10-AP waveform}$$

$$\mathcal{L} = \mathcal{L}_\text{half} + \mathcal{L}_\text{1ap} + \mathcal{L}_\text{10ap}$$

Waveforms are normalized to their peak before computing MSE so that 1 AP (0.28 ΔF/F) and 10 AP (4.7 ΔF/F) contribute equally. The half-time term uses log-scale because the empirical values span 0.5–3.0 s (6× range).

---

## What each term constrains

| Term | Primary parameters driven |
|------|--------------------------|
| `L_half` | `τ`, `τr`, `g` (saturation lengthens decay at high AP counts) |
| `L_1ap` | `A`, `τr`, `n`, `c0` (single-AP, unsaturated regime) |
| `L_10ap` | `g`, `n`, `c0` (calcium accumulation activates saturation and cooperativity) |

---

## Porting checklist

To use a different optimizer, provide:

1. **6 continuous parameters** within the bounds above.
2. **`forward_sim(τ, τr, A, g, c0, n)`** — runs `calcium_trace` + `delta_f_over_f` for
   `[1, 2, 3, 5, 10, 20, 30, 50, 100, 1, 10]` spike counts in one batched call.
3. **`compute_losses(ΔFs, t)`** — returns `L_half`, `L_1ap`, `L_10ap` as above.
4. **Scalar return** `L_total` to minimize.

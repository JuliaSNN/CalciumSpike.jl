# Equivalence check between the MATLAB deconvolution pipeline used for the
# real imaging data in the optimization project
# (`../GCaMP6s/CalciumImagingPipeline_Deconvolve.m`) and CalciumSpike.jl's own
# deconvolution (`calcium_postprocess`).
#
# Not part of the CalciumSpike module — a standalone verification script:
# 1. Generate a synthetic spike train.
# 2. Run it through CalciumSpike's forward model to get a synthetic ΔF/F-like
#    fluorescence trace (the same kind of signal the MATLAB pipeline expects
#    from suite2p).
# 3. Deconvolve that SAME trace two independent ways: the MATLAB pipeline
#    (called live from Julia via MATLAB.jl) and CalciumSpike's own
#    `calcium_postprocess`.
# 4. Compare the two deconvolved traces (correlation + overlay plot).
#
# High correlation between the two outputs is the equivalence criterion —
# the two pipelines use different baseline/smoothing implementations, so
# exact numeric agreement isn't expected or required.

using CalciumSpike
using SNNModels
using MATLAB
using CairoMakie
using Statistics
@load_units

CalciumSpike.activate_MLSpike()
mat_dir = joinpath(@__DIR__, "..", "GCaMP6s")
mat"addpath($mat_dir)"

# ── Synthetic ground truth ──────────────────────────────────────────────────
# GCaMP6s forward-model parameters (Chen et al. 2013, Nature 499:295–300),
# taken verbatim from the real production config used throughout the
# SpatioTemporalTransformation paper:
# `papers/SpatioTemporalTransformation/library/parameters/params_analysis.jl`
# (`analysis.calcium_param.forward`). σ is set to 0 here (production uses
# 0.01) to isolate the deconvolution step from noise-estimation effects —
# everything else matches production exactly.
params = CaModel(
    τ  = 1.83s,
    τr = 197ms,
    A  = 0.3,
    g  = 0.05,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.0,   # production: 0.01 — zeroed here to isolate deconvolution
    c0 = 0.23,
    n  = 2.05,
)

sr       = 50Hz    # matches params_analysis.jl's calcium_param.sr exactly
sim_time = 60s

spikes = generate_spike_train(ntrial = 1, T = sim_time / s, rate = 2.0)
gcamp  = calcium_trace(spikes[1], sr, (0s, sim_time); params)

# MATLAB expects suite2p-style F/Fneu: neurons × frames. CalciumSpike's
# forward model has no neuropil contamination, so Fneu = 0 and neuropil
# subtraction is skipped on the MATLAB side (subtract_Neuropil = false).
F    = reshape(Float64.(gcamp.F), 1, :)
Fneu = zeros(size(F))

# `tau` here is the DECONVOLUTION time constant (production's
# `calcium_param.post.τ = 2s`), not the forward model's own τ (1.83s) above —
# params_analysis.jl deliberately decouples the two.
MATLAB_POST_TAU = 2s
dt_s  = (1 / sr) / s
tau_s = MATLAB_POST_TAU / s

# ── MATLAB pipeline ──────────────────────────────────────────────────────────
FDeconv_matlab = nothing
n_frames = size(F, 2)
mat"""
[$FDeconv_matlab, ~] = CalciumImagingPipeline_Deconvolve($F, $Fneu, 1, $n_frames, $dt_s, false, $tau_s);
"""
FDeconv_matlab = vec(FDeconv_matlab)

# ── CalciumSpike's own deconvolution ────────────────────────────────────────
# `match_matlab_scale` (false reproduces CalciumSpike's own defaults,
# unmatched against MATLAB) applies every correction needed to put both
# pipelines on the same absolute scale, not just the same timing.
#
# Parameter correspondence — every value below is the real production value
# from `params_analysis.jl`'s `analysis.calcium_param.post` (a
# `CaPostProcess`), cross-checked against its MATLAB-side origin in
# CalciumImagingPipeline_Deconvolve.m / gaussianfilter.m — nothing here is
# fit or guessed:
#
# | Parameter           | MATLAB source                                     | params_analysis.jl / Julia argument       |
# |----------------------|----------------------------------------------------|---------------------------------------------|
# | Baseline window      | `FiltWindow=100` samples (fn_filt) at 31.5Hz       | `post.baseline_window = 3175ms` (100/31.5Hz) |
# | Baseline percentile  | `f0Prctle=20`                                       | `delta_f_over_f(q=0.20)` (already CalciumSpike's own default) |
# | Decay constant τ     | `tau` argument (real seconds)                       | `post.τ = 2s` (decoupled from the forward model's own τ=1.83s above) |
# | Smoothing sigma      | `gaussianfilter(temp,1,1.5,1)` — 1.5 samples at 31.5Hz | `post.σsmooth = 47.6ms` (1.5/31.5Hz)   |
# | Smoothing causality  | `gaussianfilter(temp,1,1.5,1)` — 4th arg (forward)  | `post.skewed = :left`                     |
# | Post-processing amplitude | not applicable (MATLAB never divides by `A`) | `post.A = 1.0` (fixed — NOT the forward model's A=0.3) |
# | Output scale         | never divided by `dt` or `A` (derivation below)     | `dec_julia .*= (MATLAB_POST_A * δT)`      |
#
# Output-scale derivation: MATLAB's `temp[i] = (V[i+1]-V[i]) + V[i]*dt/tau`
# is `dt_real * firing[i]` (never divided by `dt` or by the per-spike
# amplitude `A`). CalciumSpike's `deconvolve_df_f` computes
# `(dF/dt + F/τ)/A = firing[i]/A` — a properly normalized rate. Since `dt`/`τ`
# here are SNNModels-internal ms-tick numbers (1 s = 1000, see
# `src/SNNModels.jl/src/utils/unit.jl`) rather than real seconds, both terms
# are divided by values 1000× too large, so CalciumSpike's computed rate
# comes out exactly 1000× smaller than the real-seconds formula predicts.
# Net: both factors combine into a single `1/(A*δT)`, δT being the ms-tick
# sample spacing already in scope — so multiplying `dec_julia` by `A * δT`
# converts it into MATLAB's implicit units exactly, with no approximation.
match_matlab_scale = true

# params_analysis.jl: analysis.calcium_param.post = CaPostProcess(...)
MATLAB_BASELINE_WINDOW_MS  = 3175.0f0  # "100 frames at 31.5 Hz — Bagur et al. 2025"
MATLAB_BASELINE_PERCENTILE = 0.20f0    # CalciumImagingPipeline_Deconvolve.m: f0Prctle
MATLAB_SMOOTH_SIGMA_MS     = 47.6f0    # "1.5 frames at 31.5 Hz — Bagur et al. 2025"
MATLAB_SMOOTH_CAUSAL       = true      # gaussianfilter(temp,1,1.5,1) — 4th argument
MATLAB_POST_A              = 1.0f0     # post.A — fixed, decoupled from the forward model's A

δT = gcamp.t[2] - gcamp.t[1]   # sample spacing, ms-tick units (== 1/sr = 1/50Hz here)

baseline_window = match_matlab_scale ? MATLAB_BASELINE_WINDOW_MS : 3000.0f0
ΔF, r   = delta_f_over_f(gcamp.t, gcamp.F; baseline_window, q = MATLAB_BASELINE_PERCENTILE)
σsmooth = match_matlab_scale ? MATLAB_SMOOTH_SIGMA_MS : 100.0f0
skewed  = match_matlab_scale ? (MATLAB_SMOOTH_CAUSAL ? :left : :none) : :none
post_A  = match_matlab_scale ? MATLAB_POST_A : params.A
post      = CaPostProcess(τ = MATLAB_POST_TAU, A = post_A, σsmooth = σsmooth, skewed = skewed)
dec_julia = calcium_postprocess(ΔF, r, post)
if match_matlab_scale
    dec_julia = dec_julia .* (post_A * δT)
end

# ── Pre-fix (buggy) gaussian_smooth, for comparison only ────────────────────
# Verbatim copy of SNNModels.jl's gaussian_smooth as it was before the fix
# (commit 614cfca) — offsets compared directly against sigma in xs' physical
# units, instead of being converted to those units first. Kept local to this
# script (not reverting the real, fixed function) purely to show the
# difference the bug made.
function gaussian_smooth_buggy(xs::RT, x::T, σ::R; skewed::Symbol=:none) where {T<:AbstractVector, R<:Real, RT<:AbstractVector}
    iszero(σ) && return copy(x)
    step_x = xs[2] - xs[1]
    half = ceil(Int, 3σ / step_x)
    offsets = if skewed === :left
        -half:0
    elseif skewed === :right
        0:half
    else
        -half:half
    end
    kernel = exp.(.-Float64.(offsets).^2 ./ (2σ^2))   # bug: offsets in samples, σ in xs' units
    kernel ./= sum(kernel)
    n = length(x)
    out = similar(x)
    lo = first(offsets)
    @inbounds for i in 1:n
        s = 0.0; w = 0.0
        @fastmath for (j, k) in enumerate(kernel)
            idx = i + lo + (j - 1)
            (idx < 1 || idx > n) && continue
            s += k * x[idx]
            w += k
        end
        out[i] = s / w
    end
    return out
end

dec_exp_prebuggy = deconvolve_df_f(ΔF, 1 / (r[2] - r[1]), post.τ, post.A)
dec_julia_buggy  = gaussian_smooth_buggy(r, dec_exp_prebuggy, σsmooth; skewed)
if match_matlab_scale
    dec_julia_buggy = dec_julia_buggy .* (post_A * δT)
end

# ── Compare ──────────────────────────────────────────────────────────────────
n = min(length(FDeconv_matlab), length(dec_julia), length(dec_julia_buggy))
c       = cor(FDeconv_matlab[1:n], dec_julia[1:n])
c_buggy = cor(FDeconv_matlab[1:n], dec_julia_buggy[1:n])
@info "Correlation MATLAB vs CalciumSpike (fixed gaussian_smooth): $(round(c, digits = 3))"
@info "Correlation MATLAB vs CalciumSpike (pre-fix, buggy gaussian_smooth): $(round(c_buggy, digits = 3))"

t_s = collect(r[1:n]) ./ 1000

fig = Figure(size = (900, 750))
windows = [(0, 20), (20, 40), (40, 60)]
for (row, (t0, t1)) in enumerate(windows)
    ax = Axis(fig[row, 1];
              xlabel = row == length(windows) ? "Time (s)" : "",
              ylabel = "Deconvolved activity (a.u.)",
              title  = row == 1 ?
                       "MATLAB vs CalciumSpike.jl deconvolution (fixed r=$(round(c, digits=3)), buggy r=$(round(c_buggy, digits=3)))" : "")
    lines!(ax, t_s, FDeconv_matlab[1:n];    label = "MATLAB (CalciumImagingPipeline_Deconvolve)")
    lines!(ax, t_s, dec_julia[1:n];         label = "CalciumSpike.jl (fixed gaussian_smooth)")
    lines!(ax, t_s, dec_julia_buggy[1:n];   label = "CalciumSpike.jl (pre-fix, buggy gaussian_smooth)", linestyle = :dash)
    xlims!(ax, t0, t1)
    row == 1 && axislegend(ax)
end

save(joinpath(@__DIR__, "matlab_vs_calciumspike.png"), fig)
fig

##
fig = Figure(size = (900, 750))
windows = [(0, 20), (45, 48), (40, 60)]
for (row, (t0, t1)) in enumerate(windows)
    ax = Axis(fig[row, 1];
              xlabel = row == length(windows) ? "Time (s)" : "",
              ylabel = "Deconvolved activity (a.u.)",
              title  = row == 1 ?
                       "MATLAB vs CalciumSpike.jl deconvolution (fixed r=$(round(c, digits=3)), buggy r=$(round(c_buggy, digits=3)))" : "")
    lines!(ax, t_s, FDeconv_matlab[1:n];    label = "MATLAB (CalciumImagingPipeline_Deconvolve)")
    lines!(ax, t_s, dec_julia[1:n];         label = "CalciumSpike.jl (fixed gaussian_smooth)")
    lines!(ax, t_s, dec_julia_buggy[1:n];   label = "CalciumSpike.jl (pre-fix, buggy gaussian_smooth)", linestyle = :dash)
    xlims!(ax, t0, t1)
    row == 1 && axislegend(ax)
end
fig


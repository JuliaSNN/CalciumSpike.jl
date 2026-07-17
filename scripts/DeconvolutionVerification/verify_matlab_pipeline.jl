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
params = CaModel(
    τ  = 2s,
    A  = 0.2,
    F0 = 1.0,
    η  = 0.0,
    σ  = 0.0,   # noise-free: isolate the deconvolution step itself
)

sr       = 50Hz
sim_time = 60s

spikes = generate_spike_train(ntrial = 1, T = sim_time / s, rate = 2.0)
gcamp  = calcium_trace(spikes[1], sr, (0s, sim_time); params)

# MATLAB expects suite2p-style F/Fneu: neurons × frames. CalciumSpike's
# forward model has no neuropil contamination, so Fneu = 0 and neuropil
# subtraction is skipped on the MATLAB side (subtract_Neuropil = false).
F    = reshape(Float64.(gcamp.F), 1, :)
Fneu = zeros(size(F))

dt_s  = (1 / sr) / s
tau_s = params.τ / s

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
# Parameter correspondence — every MATLAB-side value below is read directly
# from CalciumImagingPipeline_Deconvolve.m / gaussianfilter.m, not fit or
# guessed:
#
# | Parameter           | MATLAB source                                     | Julia argument                          |
# |----------------------|----------------------------------------------------|------------------------------------------|
# | Baseline window      | `FiltWindow=100` (samples, fn_filt)                | `delta_f_over_f(baseline_window=100δT)`   |
# | Baseline percentile  | `f0Prctle=20`                                       | `delta_f_over_f(q=0.20)` (already CalciumSpike's own default) |
# | Decay constant τ     | `tau` argument (real seconds)                       | `CaPostProcess(τ=params.τ)`               |
# | Smoothing sigma      | `gaussianfilter(temp,1,1.5,1)` — 3rd arg (samples)  | `CaPostProcess(σsmooth=1.5δT)`            |
# | Smoothing causality  | `gaussianfilter(temp,1,1.5,1)` — 4th arg (forward)  | `CaPostProcess(skewed=:left)`             |
# | Output scale         | never divided by `dt` or `A` (derivation below)     | `dec_julia .*= (A * δT)`                  |
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

MATLAB_FILT_WINDOW_SAMPLES  = 100     # CalciumImagingPipeline_Deconvolve.m: FiltWindow
MATLAB_BASELINE_PERCENTILE  = 0.20f0  # CalciumImagingPipeline_Deconvolve.m: f0Prctle
MATLAB_SMOOTH_SIGMA_SAMPLES = 1.5     # gaussianfilter(temp,1,1.5,1) — 3rd argument
MATLAB_SMOOTH_CAUSAL        = true    # gaussianfilter(temp,1,1.5,1) — 4th argument

δT = gcamp.t[2] - gcamp.t[1]   # sample spacing, ms-tick units (== 1/sr)

baseline_window = match_matlab_scale ? MATLAB_FILT_WINDOW_SAMPLES * δT : 3000.0f0
ΔF, r   = delta_f_over_f(gcamp.t, gcamp.F; baseline_window, q = MATLAB_BASELINE_PERCENTILE)
σsmooth = match_matlab_scale ? MATLAB_SMOOTH_SIGMA_SAMPLES * δT : 100.0f0
skewed  = match_matlab_scale ? (MATLAB_SMOOTH_CAUSAL ? :left : :none) : :none
post      = CaPostProcess(τ = params.τ, A = params.A, σsmooth = σsmooth, skewed = skewed)
dec_julia = calcium_postprocess(ΔF, r, post)
if match_matlab_scale
    dec_julia = dec_julia .* (params.A * δT)
end

# ── Compare ──────────────────────────────────────────────────────────────────
n = min(length(FDeconv_matlab), length(dec_julia))
c = cor(FDeconv_matlab[1:n], dec_julia[1:n])
@info "Correlation between MATLAB and CalciumSpike deconvolution: $(round(c, digits = 3))"

t_s = collect(r[1:n]) ./ 1000

fig = Figure(size = (900, 750))
windows = [(0, 20), (20, 40), (40, 60)]
for (row, (t0, t1)) in enumerate(windows)
    ax = Axis(fig[row, 1];
              xlabel = row == length(windows) ? "Time (s)" : "",
              ylabel = "Deconvolved activity (a.u.)",
              title  = row == 1 ?
                       "MATLAB vs CalciumSpike.jl deconvolution (r = $(round(c, digits = 3))" *
                       (match_matlab_scale ? ", scale-matched)" : ", raw)") : "")
    lines!(ax, t_s, FDeconv_matlab[1:n]; label = "MATLAB (CalciumImagingPipeline_Deconvolve)")
    lines!(ax, t_s, dec_julia[1:n];      label = "CalciumSpike.jl (calcium_postprocess)")
    xlims!(ax, t0, t1)
    row == 1 && axislegend(ax)
end

save(joinpath(@__DIR__, "matlab_vs_calciumspike.png"), fig)
fig


# debug_sigma_leakage.jl
#
# Show how symmetric Gaussian smoothing after calcium deconvolution leaks
# information across epoch boundaries in a 200-neuron synthetic population.
#
# Firing rate schedule (heterogeneous Pareto-distributed across neurons):
#   0 – 500 ms :    5 Hz mean  (baseline)
#   500 – 650 ms :  15 Hz mean  (onset burst)
#   650 – 1000 ms :  5 Hz mean  (sustained baseline)
#   1000 – 1150 ms: 10 Hz mean  (offset burst)
#   1150 – 1500 ms : 5 Hz mean  (baseline)
#
# Pipeline: spikes → calcium_trace (CaModel) → ΔF/F → deconvolve → gaussian_smooth(σ)
# σ sweep: [1, 2, 5, 10, 50, 100, 200, 500] ms
#
# Note: internal time unit = ms  (s = 1000, ms = 1, Hz = 0.001)
#
# Output: scripts/GCaMP6s/debug_sigma_leakage.png

using CalciumSpike
using SNNModels
@load_units
using CairoMakie, Statistics, Random
using Distributions: Pareto, Poisson

# ── Calcium parameters (from STT params_analysis) ─────────────────────────────
const CA_MODEL = CaModel(
    τ  = 1.83f0s,
    τr = 197f0ms,
    A  = 0.3f0,
    g  = 0.05f0,
    F0 = 1.0f0,
    η  = 0.0f0,
    σ  = 0.01f0,
    c0 = 0.23f0,
    n  = 2.05f0,
)

const CA_POST = CaPostProcess(
    τ       = 2f0s,
    A       = 0.2f0,
    σsmooth = 100f0ms,   # default; swept below
    skewed  = :right,
)

const SR         = 50Hz           # imaging sampling rate (0.05 samples/ms → dt=20ms)
const STAB_MS    = 5000f0         # 5 s stabilization before stimulus
const STIM_MS    = 1500f0         # stimulus window length
const INTERVAL   = (0f0, STAB_MS + STIM_MS)
const N_NEURONS  = 200

# ── Firing rate epochs (t_start_ms, t_end_ms, mean_rate_Hz) ───────────────────
# First epoch is the 5-second stabilization at baseline rate.
# Stimulus epochs are offset by STAB_MS.
const EPOCHS = [
    (0f0,                 STAB_MS,                5f0),   # 5 s stabilization
    (STAB_MS + 0f0,    STAB_MS + 500f0,   5f0),
    (STAB_MS + 500f0,  STAB_MS + 650f0,  15f0),
    (STAB_MS + 650f0,  STAB_MS + 1000f0,  5f0),
    (STAB_MS + 1000f0, STAB_MS + 1150f0, 10f0),
    (STAB_MS + 1150f0, STAB_MS + 1500f0,  5f0),
]

# ── Generate synthetic population spikes ───────────────────────────────────────
# Each neuron has a Pareto rate multiplier (α=2, θ=0.5 → mean=1).
# Expected spikes per epoch: rate_Hz * duration_s (Poisson draw).
rng = MersenneTwister(42)
rate_mult = rand(rng, Pareto(2f0, 0.5f0), N_NEURONS)   # mean multiplier = 1.0

function generate_spikes(n_neurons, epochs, rate_mult, rng)
    spikes = [Float32[] for _ in 1:n_neurons]
    for (t0_ms, t1_ms, r_hz) in epochs
        dur_ms = t1_ms - t0_ms
        dur_s  = dur_ms / s                    # convert ms → seconds for rate calc
        for n in 1:n_neurons
            rate_hz = r_hz * rate_mult[n]
            n_sp    = rand(rng, Poisson(rate_hz * dur_s))
            append!(spikes[n], Float32.(t0_ms .+ dur_ms .* sort(rand(rng, n_sp))))
        end
    end
    return sort!.(spikes)
end

spikes = generate_spikes(N_NEURONS, EPOCHS, rate_mult, rng)
@info "Spikes generated" n_neurons=N_NEURONS total_spikes=sum(length.(spikes))

# ── Forward calcium model ──────────────────────────────────────────────────────
@info "Running calcium forward model..."
Fs, t = calcium_trace(spikes, SR, INTERVAL; params = CA_MODEL)
@info "calcium_trace done" n_neurons=length(Fs) n_timepoints=length(t) t_start=t[1] t_end=t[end]
# Use the stabilization period to estimate F0 baseline
ΔFs, t = delta_f_over_f(t, Fs)
t_ms = collect(t)   # already in ms

# Indices for the stimulus window only (trim stabilization for display)
stim_idx = findall(>=(STAB_MS), t_ms)
t_plot   = t_ms[stim_idx] .- STAB_MS   # relative time: 0 = stimulus onset

# ── Sigma sweep ────────────────────────────────────────────────────────────────
sigmas_ms = [1, 2, 5, 10, 50, 100, 200, 500]
@info "Sweeping σsmooth: $(sigmas_ms) ms"

pop_mean = map(sigmas_ms) do σ_ms
    @info "  σ = $(σ_ms) ms"
    post = CaPostProcess(τ = CA_POST.τ, A = CA_POST.A, σsmooth = Float32(σ_ms) * ms)
    dec  = calcium_postprocess(ΔFs, t_ms, post)  # Vector{Vector{Float32}}
    mat  = reduce(hcat, dec)                      # (n_timepoints × n_neurons)
    full = vec(mean(mat, dims = 2))
    full[stim_idx]                                # trim to stimulus window
end

# ── Figure ─────────────────────────────────────────────────────────────────────
ncols   = 4
nrows   = 2
palette = cgrad(:roma, length(sigmas_ms), categorical = true)

fig = Figure(size = (ncols * 340, nrows * 270))

epoch_bounds = [0, 500, 650, 1000, 1150, 1500]
epoch_colors = [:grey80, :tomato, :grey80, :steelblue, :grey80]

for (idx, (σ, trace)) in enumerate(zip(sigmas_ms, pop_mean))
    row = (idx - 1) ÷ ncols + 1
    col = (idx - 1) % ncols + 1
    ax  = Axis(fig[row, col];
        title  = "σ = $(σ) ms",
        xlabel = row == nrows ? "Time (ms)" : "",
        ylabel = col == 1 ? "Deconv. activity (a.u.)" : "",
    )
    for (i, (lo, hi)) in enumerate(zip(epoch_bounds[1:end-1], epoch_bounds[2:end]))
        vspan!(ax, Float32(lo), Float32(hi); color = (epoch_colors[i], 0.18))
    end
    for tb in epoch_bounds[2:end-1]
        vlines!(ax, Float32(tb); color = :black, linestyle = :dash, linewidth = 0.8)
    end
    lines!(ax, t_plot, trace; color = palette[idx], linewidth = 1.8)
    xlims!(ax, 0, 1500)
end

Label(fig[0, :];
    text     = "Calcium σsmooth leakage — population mean (200 neurons, Pareto heterogeneous rates)",
    fontsize = 13, font = :bold)

# ── Save ───────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "debug_sigma_leakage.png")
save(out, fig; px_per_unit = 2)
@info "Saved: $out"

fig
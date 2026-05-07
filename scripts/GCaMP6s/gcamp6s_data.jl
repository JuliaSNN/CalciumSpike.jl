# GCaMP6s empirical data — Chen et al. 2013, Nature 499, 295–300
# "Ultrasensitive fluorescent proteins for imaging neuronal activity"
#
# Figure 1b — Fluorescence changes averaged across neurons and wells for GCaMP6s.
#   Top panel:    response to 1 action potential
#   Bottom panel: response to 10 action potentials
#
# Figure 1f — Half decay time as a function of stimulus strength (AP count) for GCaMP6s.
#   n = 11 wells; error bars = s.e.m.
#
# Digitized from published figures. Time axes are relative to stimulus onset (seconds).
# ΔF/F₀ amplitudes are as reported (not normalized).

# --- Figure 1f: half-decay time vs AP count ---
# Digitized from plot-data_halftime.csv. x values are digitizer output (≈ AP count).

const DATA_X = Int[1, 2, 3, 5, 10, 20, 40, 80]

# Half-decay time after stimulus (seconds) — GCaMP6s
const T_HALF_EMPIRICAL = Float32[1.029, 1.428, 1.493, 1.581, 1.809, 2.267, 2.730, 3.276]

# --- Figure 1b top: 1 AP transient ---
# Time (s) relative to spike, ΔF/F₀

const T_1AP = Float32[
    0.054, 0.114, 0.144, 0.281, 0.569, 0.737, 0.940,
    1.108, 1.263, 1.449, 1.701, 1.952, 2.198, 2.455, 2.719, 2.964,
]
const DF_1AP = Float32[
    0.105, 0.168, 0.209, 0.261, 0.253, 0.230, 0.205,
    0.190, 0.168, 0.148, 0.119, 0.096, 0.079, 0.071, 0.056, 0.041,
]

const AP_1 = 1
const AP_10 = 10
# const T_SPIKE = 500f0   # ms, time of spike(s) in simulation

# --- Figure 1b bottom: 10 AP transient ---
# Time (s) relative to stimulus offset, ΔF/F₀

const T_10AP = Float32[
    -0.249, -0.006, 0.077, 0.127, 0.215, 0.387, 0.492, 0.641,
     0.829,  1.017, 1.221, 1.403, 1.580, 1.762, 1.961, 2.127,
     2.287,  2.453, 2.646, 2.862,
]
const DF_10AP = Float32[
    -0.013,  0.013, 0.971, 2.297, 3.556, 4.383, 4.541, 4.528,
     4.409,  4.213, 3.963, 3.648, 3.412, 3.176, 2.887, 2.717,
     2.520,  2.323, 2.139, 1.929,
]
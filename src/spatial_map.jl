# ── Generic helpers ───────────────────────────────────────────────────────────

function neuron_coords(points, neurons)
    if isa(points, Matrix)
        points[neurons, 1], points[neurons, 2]
    else
        [p[1] for p in points][neurons], [p[2] for p in points][neurons]
    end
end

"""
    compute_spatial_map(activity, time_indices, neurons, points, recordings; kwargs...)

Compute spatial activity maps for given neuron activity matrix and time intervals.

# Arguments
- `activity`: neurons × time matrix
- `intervals`: vector of `(t_start, t_end)` pairs (ms)
- `neurons`: indices into `points`
- `points`: spatial positions (Matrix or Vector of tuples)
- `recordings`: used to convert ms intervals to time indices
"""
function compute_spatial_map(;activity, baseline_activity, entries, time_indices, neurons, points, kwargs...)
    pre_results  = Dict{Symbol, Array{Float32, 3}}()
    results = Dict{Symbol, Array{Float32, 3}}()
    _m = mean(baseline_activity, dims = 2)[:, 1]
    _s = std(baseline_activity,  dims = 2)[:, 1]
    xs, ys = neuron_coords(points, neurons)

    map(keys(entries) |> collect) do entry
        recording_keys = entries[entry]
        _activity = mean([activity[kk] for kk in recording_keys])
        _activity = clamp.((_activity .- _m) ./ (_s), -10, 10)
        _activity[isnan.(_activity)] .= 0.0
        pre_results[entry], _, _ = SNNModels.spatial_activity(
                (xs, ys),
                _activity[neurons, :];
                kwargs...,
                T = time_indices,
            )
    end
    # return results
    mean_act = mean([x[:,:,1:1] for x in values(pre_results)])
    for key in keys(pre_results)
        mean_act = pre_results[key][:,:,1:1] 
        results[key] = pre_results[key][:,:,2:end] .- mean_act
    end
    return results
end

"""
    spatial_map_plot!(fig, activity_dict, ordered_keys, x_range, y_range; kwargs...)

Generic contourf grid: rows = `ordered_keys`, cols = time windows in each activity array.

# Keyword arguments
- `row_labels`: label per row (default: `string.(ordered_keys)`)
- `col_labels`: header label per column at row 0 (optional)
- `x_axis_label`, `y_axis_label`, `row_group_label`: axis labels
- `last_row_xticks`: xticks tuple applied to last row axes only (optional)
- `folder`, `filename`: both required to save
"""
function spatial_map_plot!(
    fig, activity_dict, ordered_keys, x_range, y_range;
    row_labels      = string.(ordered_keys),
    col_labels      = nothing,
    x_axis_label    = "x",
    y_axis_label    = "y",
    row_group_label = nothing,
    last_row_xticks = nothing,
    folder          = nothing,
    filename        = nothing,
    kwargs...,
)
    levels = range(-1, 1, 16)
    n_keys = length(ordered_keys)
    n_cols = size(activity_dict[ordered_keys[1]], 3)

    for (k, key) in enumerate(ordered_keys)
        result = activity_dict[key]
        for n in 1:n_cols
            ax = Axis(fig[k, n]; xlabel = "", ylabel = "")
            r  = result[:, :, n]
            contourf!(ax, x_range, y_range, r; levels, colormap = :balance)
            if k == n_keys
                !isnothing(last_row_xticks) && (ax.xticks = last_row_xticks)
                hideydecorations!(ax)
            else
                hidexdecorations!(ax)
                hideydecorations!(ax)
            end
        end
        Label(fig[k, n_cols + 1], row_labels[k], rotation = -pi / 2)
        rowsize!(fig, k, Fixed(50))
    end

    if !isnothing(col_labels)
        for (n, lbl) in enumerate(col_labels)
            Label(fig[0, n], lbl)
        end
    end
    Label(fig[n_keys + 1, 1:n_cols], x_axis_label)
    Label(fig[1:n_keys, 0], y_axis_label, rotation = pi / 2)
    !isnothing(row_group_label) &&
        Label(fig[1:n_keys, n_cols + 2], row_group_label, rotation = -pi / 2)
    [colsize!(fig, f, Relative(1 / n_cols)) for f in 1:n_cols]
    !isnothing(folder) && !isnothing(filename) &&
        savefig(fig, joinpath(folder, filename))
    return fig
end

export compute_spatial_map, spatial_map_plot!, neuron_coords

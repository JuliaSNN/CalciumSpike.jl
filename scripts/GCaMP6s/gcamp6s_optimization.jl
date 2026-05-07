using DrWatson
findproject() |> quickactivate
include(projectdir("loaders", "cluster.jl"))
include(projectdir("loaders", "filesystem.jl"))

@everywhere begin
    include(joinpath(@__DIR__, "gcamp6s_data.jl"))
    include(joinpath(@__DIR__, "gcamp_response_analysis.jl"))

    root    = data_root("helix")
    subpath = joinpath("data", "models", "calcium", "gcamp6fit")
    db_name = "CaModel_TPE"

    trials_path   = joinpath(root, subpath, "trials") |> mkpath
    study_path    = joinpath(root, subpath, db_name * ".db")
    artifact_path = joinpath(root, subpath, "artifacts", db_name) |> mkpath

    study = (
        name          = "gcamp6_CaModel_fit",
        db_name       = db_name,
        storage       = study_path,
        trials_path   = trials_path,
        artifact_path = artifact_path,
        root          = root,
        path          = subpath,
    )

    storage = journal_storage(study.storage; lock = JournalFileOpenLock)

    directions   = pylist(["minimize"])
    metric_names = pylist(["loss"])
    opt_study = OptunaLoader.optuna[].create_study(
        sampler = OptunaLoader.optuna[].samplers.TPESampler(
            multivariate     = true,
            constant_liar    = true,
            n_startup_trials = 400,
            consider_prior   = true,
        ),
        directions     = directions,
        study_name     = study.name,
        storage        = storage,
        load_if_exists = true,
        pruner = OptunaLoader.optuna[].pruners.PercentilePruner(
            25.0,
            n_startup_trials = 5,
            n_warmup_steps   = 0,
            interval_steps   = 1,
        ),
    )
    opt_study.set_metric_names(metric_names)

    base_config = (; study)

    struct ObjectiveGCaMP6
        config::NamedTuple
    end

    function (obj::ObjectiveGCaMP6)(trial)
        # τ and τr suggested in ms to stay in SNNModels unit system (1 s = 1000)
        τ  = pyconvert(Float32, trial.suggest_float("tau",   500.0, 3000.0))
        τr = pyconvert(Float32, trial.suggest_float("tau_r",  10.0,  500.0))
        A  = pyconvert(Float32, trial.suggest_float("A",       0.05,   1.0))
        g  = pyconvert(Float32, trial.suggest_float("g",       0.0,    0.3))
        c0 = pyconvert(Float32, trial.suggest_float("c0",      0.0,    1.5))
        n  = pyconvert(Float32, trial.suggest_float("n",       1.0,    4.0))

        params = CaModel(
            τ  = τ,
            τr = τr,
            A  = A,
            g  = g,
            c0 = c0,
            n  = n,
            F0 = 1.0f0,
            η  = 0.0f0,
            σ  = 0.0f0,
        )

        all_counts = vcat(DATA_X, [1, 10])
        ΔFs, t = sim_traces(all_counts, params)

        n_halftime = length(DATA_X)
        sim_half = [half_decay_time(ΔFs[i], t) for i in 1:n_halftime]
        sim_half_clamped = clamp.(sim_half, 1f-4, 1f3)
        loss_half = mean((log.(sim_half_clamped) .- log.(T_HALF_EMPIRICAL)) .^ 2)

        loss_1ap  = waveform_loss(ΔFs[end-1], t, T_1AP,  DF_1AP)
        loss_10ap = waveform_loss(ΔFs[end],   t, T_10AP, DF_10AP)

        total_loss = loss_half + loss_1ap + loss_10ap

        trial.set_user_attr("loss_half",  loss_half)
        trial.set_user_attr("loss_1ap",   loss_1ap)
        trial.set_user_attr("loss_10ap",  loss_10ap)
        trial.set_user_attr("tau_s",      τ  / 1000)
        trial.set_user_attr("tau_r_s",    τr / 1000)

        return total_loss
    end
end

##
@sync @distributed for _ in 1:5
    gcamp_objective = ObjectiveGCaMP6(base_config)
    opt_study.optimize(gcamp_objective, n_trials = 500)
end

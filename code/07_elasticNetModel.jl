using CSV
using DataFrames
using GLMNet
using Statistics
using Random
using StatsBase
using Tables

# set seed
import Random
Random.seed!(66)

# data

topology = CSV.read("data/cleaned/all_networks.csv", DataFrame)

cluster_medoids = CSV.read("data/outputs/cluster_medoids.csv", DataFrame)[:,1]
pc_dominant_metrics = CSV.read("data/outputs/pc_dominant_metrics.csv", DataFrame)[:,1]
pc_scores_df = CSV.read("data/outputs/pc_scores_df.csv", DataFrame)

pc_names = names(pc_scores_df)

# Combine
topology_pc = hcat(topology, pc_scores_df)

stability_vars = ["robustness", "ρ", "control", "resilience"]

# functions

function rmse(y_true, y_pred)
    return sqrt(mean((y_true .- y_pred).^2))
end

function r2_score(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return 1 - (ss_res / ss_tot)
end

function kfold_indices(n, k, repeats)
    all_splits = []

    for r in 1:repeats
        idx = shuffle(1:n)
        folds = [idx[floor(Int, (i-1)*n/k)+1 : floor(Int, i*n/k)] for i in 1:k]
        push!(all_splits, collect(folds))
    end

    return all_splits
end

function run_stability_enet(target_var, predictor_names, data_full)

    # -------------------------
    # Prepare data
    # -------------------------
    cols = vcat([target_var], predictor_names)
    df = dropmissing(data_full[:, cols])

    y_raw = Vector(df[:, target_var])
    X_raw = Matrix(df[:, predictor_names])
    n = size(X_raw, 1)

    alphas = 0:0.25:1

    # CV tracking
    best_r2 = -Inf
    best_rmse = Inf
    best_alpha = nothing
    best_lambda_value = nothing

    folds_all = kfold_indices(n, 5, 10)

    # -------------------------
    # Cross-validation loop
    # -------------------------
    for α in alphas

        lambda_scores = Dict{Float64, Vector{Float64}}()
        lambda_rmse   = Dict{Float64, Vector{Float64}}()

        for folds in folds_all
            for test_idx in folds

                train_idx = setdiff(1:n, test_idx)

                X_train = X_raw[train_idx, :]
                y_train = y_raw[train_idx]

                X_test = X_raw[test_idx, :]
                y_test = y_raw[test_idx]

                # -------------------------
                # Standardise train only
                # -------------------------
                μx = mean(X_train, dims=1)
                σx = std(X_train, dims=1)
                σx[σx .== 0] .= 1.0

                μy = mean(y_train)
                σy = std(y_train)
                σy = σy == 0 ? 1.0 : σy

                X_train_std = (X_train .- μx) ./ σx
                X_test_std  = (X_test  .- μx) ./ σx
                y_train_std = (y_train .- μy) ./ σy
                y_test_std  = (y_test  .- μy) ./ σy

                # -------------------------
                # Fit elastic net for this alpha (full lambda path)
                # -------------------------
                fold_model = glmnet(X_train_std, y_train_std, alpha=α)
                λ_path = fold_model.lambda

                # Evaluate each lambda
                for (idx, λ_val) in enumerate(λ_path)
                    y_pred = GLMNet.predict(fold_model, X_test_std, idx)[:,1]

                    r2 = r2_score(y_test_std, y_pred)
                    e  = rmse(y_test_std, y_pred)

                    push!(get!(lambda_scores, λ_val, Float64[]), r2)
                    push!(get!(lambda_rmse, λ_val, Float64[]), e)
                end
            end
        end

        # -------------------------
        # Pick best lambda for this alpha
        # -------------------------
        for (λ_val, r2_list) in lambda_scores
            mean_r2 = mean(r2_list)
            if mean_r2 > best_r2
                best_r2 = mean_r2
                best_rmse = mean(lambda_rmse[λ_val])
                best_alpha = α
                best_lambda_value = λ_val
            end
        end
    end

    # -------------------------
    # Fit final model on full dataset
    # -------------------------
    μx = mean(X_raw, dims=1)
    σx = std(X_raw, dims=1)
    σx[σx .== 0] .= 1.0
    μy = mean(y_raw)
    σy = std(y_raw)
    σy = σy == 0 ? 1.0 : σy

    X_std = (X_raw .- μx) ./ σx
    y_std = (y_raw .- μy) ./ σy

    final_model = glmnet(X_std, y_std, alpha=best_alpha)

    # Map CV-selected lambda to closest in final model
    final_lambda_idx = findmin(abs.(final_model.lambda .- best_lambda_value))[2]

    # Extract coefficients
    betas = final_model.betas[:, final_lambda_idx]

    # -------------------------
    # Build coefficient DataFrame
    # -------------------------
    coef_df = DataFrame(
        term = predictor_names,
        estimate = betas
    )

    corrs = [cor(X_raw[:, i], y_raw) for i in 1:size(X_raw,2)]
    coef_df.correlation = corrs
    coef_df.variance_explained = coef_df.estimate .* coef_df.correlation
    coef_df.metric .= target_var

    return (
        alpha = best_alpha,
        lambda_index = final_lambda_idx,
        lambda = final_model.lambda[final_lambda_idx],
        r2 = best_r2,
        rmse = best_rmse,
        coef_df = coef_df
    )
end

rep_list = Dict(
    "medoids" => cluster_medoids,
    "dominant" => pc_dominant_metrics,
    "pca_score" => pc_names,
    "complexity" => ["complexity"],
    "connectance" => ["connectance"],
    "trophicCoherence" => ["trophicCoherence"]
)

results = DataFrame(
    metric = String[],
    rep_name = String[],
    alpha = Float64[],
    lambda = Float64[],
    r2 = Float64[],
    rmse = Float64[]
)

all_coefs = DataFrame(
    metric = String[],
    rep_name = String[],
    term = String[],
    estimate = Float64[],
    correlation = Float64[],
    variance_explained = Float64[]
)

for metric in stability_vars
    for (rep_name, cols) in rep_list

        res = run_stability_enet(metric, cols, topology_pc)

        # summary results
        push!(results, (
            metric,
            rep_name,
            res.alpha,
            res.lambda,
            res.r2,
            res.rmse
        ))

        # coefficients
        coef_df = res.coef_df
        coef_df.rep_name .= rep_name   # add group label

        append!(all_coefs, coef_df)
    end
end

CSV.write("data/outputs/elasticNet_summary.csv", results)
CSV.write("data/outputs/elasticNet_coefficients.csv", all_coefs)
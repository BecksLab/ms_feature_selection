using CSV
using DataFrames
using GLMNet
using Statistics
using Random
using StatsBase
using Tables

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
    return cor(y_true, y_pred)^2
end

function standardise(df::DataFrame)
    df_std = copy(df)
    for col in names(df)
        μ = mean(skipmissing(df[!, col]))
        σ = std(skipmissing(df[!, col]))
        df_std[!, col] = (df[!, col] .- μ) ./ σ
    end
    return df_std
end

function kfold_indices(n, k, repeats)
    all_splits = []

    for r in 1:repeats
        idx = shuffle(1:n)
        folds = Iterators.partition(idx, ceil(Int, n/k))
        push!(all_splits, collect(folds))
    end

    return all_splits
end

function run_stability_enet(target_var, predictor_names, data_full)

    cols = vcat([target_var], predictor_names)
    df = data_full[:, cols]
    df = dropmissing(df)

    # Standardise
    df_std = standardise(df)

    y = df_std[:, target_var]
    X_df = df_std[:, predictor_names]
    X = Matrix(X_df)

    n = size(X, 1)

    # Grid
    alphas = 0:0.25:1
    lambdas = 10 .^ range(-4, -1, length=30)

    best_r2 = -Inf
    best_rmse = Inf
    best_alpha = nothing
    best_lambda = nothing

    folds_all = kfold_indices(n, 5, 10)

    for α in alphas
        for λ in lambdas

            r2_scores = Float64[]
            rmse_scores = Float64[]

            for folds in folds_all
                for test_idx in folds

                    train_idx = setdiff(1:n, test_idx)

                    X_train = X[train_idx, :]
                    y_train = y[train_idx]

                    X_test = X[test_idx, :]
                    y_test = y[test_idx]

                    model = glmnet(X_train, y_train,
                                   alpha=α,
                                   lambda=[λ])

                    y_pred = GLMNet.predict(model, X_test)[:,1]

                    push!(r2_scores, r2_score(y_test, y_pred))
                    push!(rmse_scores, rmse(y_test, y_pred))
                end
            end

            mean_r2 = mean(r2_scores)

            if mean_r2 > best_r2
                best_r2 = mean_r2
                best_rmse = mean(rmse_scores)
                best_alpha = α
                best_lambda = λ
            end
        end
    end

    # Final model
    model = glmnet(X, y,
                   alpha=best_alpha,
                   lambda=[best_lambda])

    # -------------------------
    # Coefficients
    # -------------------------
    betas = model.betas[:, 1]

    coef_df = DataFrame(
        term = predictor_names,
        estimate = betas
    )

    # -------------------------
    # Variance explained
    # -------------------------
    corrs = [cor(X_df[!, col], y) for col in predictor_names]

    coef_df.correlation = corrs
    coef_df.variance_explained = coef_df.estimate .* coef_df.correlation

    # add identifiers (glow_up)
    coef_df.metric .= target_var

    return (
        alpha = best_alpha,
        lambda = best_lambda,
        r2 = best_r2,
        rmse = best_rmse,
        coef_df = coef_df
    )
end

rep_list = Dict(
    "medoids" => cluster_medoids,
    "dominant" => pc_dominant_metrics,
    "pca_score" => pc_names,
    "complexity" => ["complexity"]
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
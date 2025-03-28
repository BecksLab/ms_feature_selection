using CSV
using DataFrames
using LinearAlgebra
using MultivariateStats
using Normalization
using Statistics
using StatsBase

## import data
network_struct = DataFrame(CSV.File("data/mangal/mangal_summary.csv"))

# select only the relevant variables
var_df = select(network_struct, Not([:id]))

## check for high number of zeros

# find prop of zeros
x = sum(x -> x == 0, Matrix(var_df), dims = 1) ./ nrow(network_struct)
# x[x .<= 0.80]
# find index of those above the cut off
high = findall(x -> x >= 0.95, x)

# remove vars with prop of zeros greater than cut off
zero_vars = [names(var_df, idx.I[2]) for idx in high]

## check for high pairwise correlation

cm = cor(Matrix(var_df))

high = findall(x -> 1 > abs(x) > 0.95, cm)

corr_vars = [[names(var_df, idx.I[1]); names(var_df, idx.I[2])] for idx in high]

# pca

# need to transpose for input reasons
pca_matrix = Array(var_df)'

# scale
X = Matrix(fit(ZScoreTransform, pca_matrix; dims = 2))
StatsBase.transform(X, pca_matrix)

# PCA:
M = fit(PCA, pca_matrix; pratio = 1, maxoutdim = 4)


proj = projection(M)
names(var_df)

df_transformed = projection(M)' * (pca_matrix .- mean(M))

h = plot(df_transformed[1, :], df_transformed[2, :], seriestype = :scatter, label = "")
plot!(xlabel = "PC1", ylabel = "PC2", framestyle = :box)
for i = 1:4
    plot!(
        [0, proj[i, 1]],
        [0, proj[i, 2]],
        arrow = true,
        label = names(var_df)[i],
        legend = :bottomleft,
    )
end
display(h)

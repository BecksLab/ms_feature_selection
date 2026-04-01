using CSV
using DataFrames
using SpeciesInteractionNetworks
using ProgressMeter
using JLD2

include("lib/_internals.jl");

# pull the network stats so we can get S and L
mangal_networks = load_object("data/mangal/mangal_networks.jlds")
wol_networks = load_object("data/weboflife/weboflife_networks.jlds")
vermaat_networks = load_object("data/vermaat_2009/networks.jlds")

networks_all = vcat(mangal_networks, wol_networks, vermaat_networks)

# df for network stats
random_topology = DataFrame(
    id = Any[],
    model = Any[],
    richness = Any[],
    links = Any[],
    connectance = Any[],
    diameter = Any[],
    complexity = Any[],
    distance = Any[],
    basal = Any[],
    top = Any[],
    intermediate = Any[],
    predpreyRatio = Any[],
    herbivory = Any[],
    omnivory = Any[],
    cannibal = Any[],
    l_S = Any[],
    GenSD = Any[],
    VulSD = Any[],
    TL = Any[],
    ChLen = Any[],
    ChSD = Any[],
    ChNum = Any[],
    path = Any[],
    LinkSD = Any[],
    S1 = Any[],
    S2 = Any[],
    S4 = Any[],
    S5 = Any[],
    ρ = Any[],
    centrality = Any[],
    loops = Any[],
    resilience = Any[],
    robustness = Any[],
    intervals = Any[],
    MaxSim = Any[],
    Clust = Any[],
    trophicCoherence = Any[],
    trophicVar = Any[],
    control = Any[],
);

# null models
# ============================================================
# Null model generators
# ============================================================

# 1. Equiprobable (connectance) null
function null_connectance(A)
    S = size(A, 1)
    L = sum(A)
    C = L / (S^2)

    A_rand = rand(S, S) .< C

    # remove self-loops
    for i in 1:S
        A_rand[i, i] = false
    end

    return Binary(A_rand)
end


# 2. Degree-based (Bascompte)
function null_degree(A)
    S = size(A, 1)
    L = sum(A)

    kout = sum(A, dims=2)  # out-degree
    kin  = sum(A, dims=1)  # in-degree

    P = (kout * kin) ./ L

    # ensure probabilities are ≤ 1
    P = min.(P, 1.0)

    A_rand = rand(S, S) .< P

    # remove self-loops
    for i in 1:S
        A_rand[i, i] = false
    end

    return Binary(A_rand)
end

# build network and get summary

@showprogress for i in 1:nrow(networks_all)

    N = render(Binary, networks_all.network[i])

    A = _get_matrix(N)
    edges = Binary(A)
    nodes = Unipartite(Symbol.(species(N)))
    N = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)

    # 1. CONNECTANCE NULL
    edges_con = null_connectance(A)
    N_con = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges_con)

    d1 = _network_summary(N_con)
    d1[:id] = networks_all.id[i]
    d1[:model] = "connectance_null"

    push!(random_topology, d1)


    # 2. DEGREE NULL
    edges_deg = null_degree(A)
    N_deg = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges_deg)

    d2 = _network_summary(N_deg)
    d2[:id] = networks_all.id[i]
    d2[:model] = "degree_null"

    push!(random_topology, d2)

end

## Write files
CSV.write("data/cleaned/randomNetworks.csv", random_topology)
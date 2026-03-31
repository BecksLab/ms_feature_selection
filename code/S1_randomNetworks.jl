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

# build network and get summary

@showprogress for i in 1:nrow(networks_all)

    N = render(Binary, networks_all.network[i])

    A = _get_matrix(N)
    edges = Binary(A)
    nodes = Unipartite(Symbol.(species(N)))
    N = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)

    m = sum(A)
    nswaps = 50 * m

    # 1. RANDOM SWAPS
    N_ran = deepcopy(N)
    for s in 1:nswaps
        swap!(N_ran)
    end

    d = _network_summary(N_ran)
    d[:id] = networks_all.id[i]
    d[:model] = "random"

    push!(random_topology, d)

    # 2. CONNECTANCE-CONSTRAINED SWAPS
    N_con = deepcopy(N)
    for s in 1:nswaps
        swap!(N_con, Connectance)
    end

    d2 = _network_summary(N_con)
    d2[:id] = networks_all.id[i]
    d2[:model] = "connectance_constrained"

    push!(random_topology, d2)

end

## Write files
CSV.write("data/cleaned/randomNetworks.csv", random_topology)
using CSV
using DataFrames
using SpeciesInteractionNetworks
using ProgressMeter
using JLD2

include("lib/_internals.jl");

# set seed
import Random
Random.seed!(66)

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

# Nice networks

# get summary csv
net_summs = CSV.read("data/cleaned/all_networks.csv", DataFrame)

n_reps = 10

niche_topology = DataFrame(
    id = Int[],
    run = Int[],
    richness = Int[],
    links = Int[],
    connectance = Float64[],
    diameter = Float64[],
    distance = Float64[],
    basal = Float64[],
    top = Float64[],
    intermediate = Float64[],
    predpreyRatio = Float64[],
    herbivory = Float64[],
    omnivory = Float64[],
    cannibal = Float64[],
    l_S = Float64[],
    GenSD = Float64[],
    VulSD = Float64[],
    TL = Float64[],
    ChLen = Float64[],
    ChSD = Float64[],
    ChNum = Float64[],
    path = Float64[],
    LinkSD = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
    centrality = Float64[],
    loops = Float64[],
    intervals = Float64[],
    MaxSim = Float64[],
    Clust = Float64[],
    trophicCoherence = Float64[],
    trophicVar = Float64[],
);

for j in 1:nrow(net_summs)

    # get network stats
    S = net_summs.richness[j]
    C = net_summs.connectance[j]

    for i in 1:n_reps

        # Build niche web
        N = structuralmodel(NicheModel, S, C)

        # metrics
        A = _get_matrix(N)
        L = links(N)
        S = richness(N)

        _gen = SpeciesInteractionNetworks.generality(N)
        gen = collect(values(_gen))
        _vul = SpeciesInteractionNetworks.vulnerability(N)
        vul = collect(values(_vul))
        ind_maxgen = findmax(gen)[2]
        l_s = L / S
        top = sum(vec(sum(A, dims = 1) .== 0))
        basal = sum(vec(sum(A, dims = 2) .== 0))
        int = (S - (basal + top))
        tl = trophic_level(N)

        chain = chain_metrics(N; max_depth=6)

        d2 = Dict{Symbol,Any}(
            :id => j,
            :run => i,
            :richness => S,
            :links => L,
            :connectance => connectance(N),
            :diameter => diameter(A),
            :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
            :basal => basal / S,
            :top => top / S,
            :intermediate => int / S,
            :herbivory => length(herbivore(N)) / S,
            :omnivory => length(omnivore(N)) / S,
            :cannibal => length(cannibal(N)) / S,
            :predpreyRatio => (basal + int)/(top + int),
            :l_S => l_s,
            :GenSD => std(gen) / l_s,
            :VulSD => std(vul) / l_s,
            :TL => mean(Float64.(values(tl))),
            :ChLen => chain.ChLen,
            :ChSD  => chain.ChSD,
            :ChNum => chain.ChNum,
            :path => mean(Float64.(pathlengths(N))),
            :LinkSD => std(values(SpeciesInteractionNetworks.degree(N))) / l_s,
            :S1 => length(findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)))/(richness(N)^2),
            :S2 => length(findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)))/(richness(N)^2),
            :S4 => length(findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)))/(richness(N)^2),
            :S5 => length(findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)))/(richness(N)^2),
            :centrality => mean(Float64.(values(centrality(N)))),
            :loops => length(loops(N)) / S,
            :intervals => intervality(A),
            :MaxSim => max_sim(N),
            :Clust => clustering(A),
            :trophicCoherence => trophic_coherence(N),
            :trophicVar => trophic_variance(N),
        )

        push!(niche_topology, d2)
    
    end
    
end

CSV.write("data/cleaned/nicheNetworks.csv", niche_topology)

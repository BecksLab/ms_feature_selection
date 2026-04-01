# General sundry internal functions

using Extinctions
using Graphs
using GraphsMatching
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics

# import other scripts with functions
include("diameter.jl")
include("intervality.jl")

"""
_network_summary(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    Returns the 'summary statistics' for a network
"""
function _network_summary(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

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

    D = Dict{Symbol,Any}(
        :richness => S,
        :links => L,
        :connectance => connectance(N),
        :diameter => diameter(A),
        :complexity => complexity(N),
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
        :TL => mean(collect(values(tl))),
        :ChLen => chain.ChLen,
        :ChSD  => chain.ChSD,
        :ChNum => chain.ChNum,
        :path => mean(pathlengths(N)),
        :LinkSD => std(values(SpeciesInteractionNetworks.degree(N))) / l_s,
        :S1 =>
            length(
                findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)),
            )/(richness(N)^2),
        :S2 =>
            length(
                findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)),
            )/(richness(N)^2),
        :S4 =>
            length(
                findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)),
            )/(richness(N)^2),
        :S5 =>
            length(
                findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)),
            )/(richness(N)^2),
        :ρ => spectralradius(N),
        :centrality => mean(collect(values(centrality(N)))),
        :loops => length(loops(N)) / S,
        :resilience => resilience(extinction(N)),
        :robustness => robustness(N; threshold = 50,
                                  remove_disconnected = true),
        :intervals => intervality(A),
        :MaxSim => max_sim(N),
        :Clust => clustering(A),
        :trophicCoherence => trophic_coherence(N),
        :trophicVar => trophic_variance(N),
        :control => structural_controllability(N)
)
    return D
end

"""
_get_matrix(N::SpeciesInteractionNetwork)

    Internal function to return a matrix of interactions from a
    SpeciesInteractionNetwork
"""
function _get_matrix(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    species = richness(N)
    n = zeros(Bool, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i, j] == true
                n[i, j] = true
            end
        end
    end

    return n
end

"""
diameter(N::SpeciesInteractionNetwork)

    Calculates the diameter of a food web. Where diameter is the longest 
    shortest path between two nodes
"""
function diameter(N::SpeciesInteractionNetwork)

    # extract species names
    spp = species(N)
    # empty vector for storing shortest path for each spp
    shortpath = zeros(Int64, length(spp))

    # get shortest path
    for i in eachindex(spp)

        paths = collect(values(shortestpath(N, spp[i])))

        if length(paths) > 0
            shortpath[i] = findmax(paths)[1]
        end
    end

    # return max shortest path
    return findmax(shortpath)[1]
end

_parser(x) = parse(Int, x)


"""
trophic_level(N::SpeciesInteractionNetwork)

    Calculates the trophic level of all species in a network using the average 
    shortest path from the prey of species 𝑖 to a basal species (prey-averaged)


"""
function trophic_level(N::SpeciesInteractionNetwork)

    A = _get_matrix(N) # Ensure A is dense for inversion.
    S = size(A, 1) # Species richness.
    out_degree = sum(A; dims = 2)
    D = -(A ./ out_degree) # Diet matrix.
    D[isnan.(D)] .= 0.0
    D[diagind(D)] .= 1.0 .- D[diagind(D)]
    # Solve with the inverse matrix.
    inverse = iszero(det(D)) ? pinv : inv
    tls = inverse(D) * ones(S)

    # create dictionary
    Dict(zip(species(N),tls))

end

"""
pathlengths(N::SpeciesInteractionNetwork)

    Returns the shortest pathlengths between all species pairs for a network
"""
function pathlengths(N::SpeciesInteractionNetwork)

    sp = species(N)
    path = Any[]

    for i in eachindex(sp)

        _path = collect(values(shortestpath(N, sp[i])))

        if length(_path) > 0
            append!(path, _path)
        end
    end

    return path
end

"""
herbivore(N::SpeciesInteractionNetwork)

    Returns a vector of species that are herbivores (only consume basal species)
"""
function herbivore(N::SpeciesInteractionNetwork)

    # find basal species
    gen = SpeciesInteractionNetworks.generality(N)
    basal = collect(keys(filter(((k, v),) -> v == 0, gen)))

    sp = species(N)
    herbivores = Any[]

    for i in eachindex(sp)

        prey = collect(successors(N, sp[i]))

        # is the prey a subset of or equal to
        if length(prey) > 0 && prey ⊆ basal
            push!(herbivores, sp[i])
        end
    end

    return herbivores
end

"""
omnivore(N::SpeciesInteractionNetwork)

    Returns a vector of species that are omnivores (feed on  > 1 species of different 
    trophic levels)
"""
function omnivore(N::SpeciesInteractionNetwork)

    omni = Any[]

    tl = trophic_level(N)
    sp = species(N)

    for i in eachindex(sp)
        prey = collect(successors(N, sp[i]))

        # return trophic level of prey
        _tls = [v for (k, v) in tl if k ∈ prey]

        if length(prey) > 1 && !allequal(_tls)
            push!(omni, sp[i])
        end
    end

    return omni
end

"""
cannibal(N::SpeciesInteractionNetwork)

    Returns a vector of species that are cannibals (feed on themselves)
"""
function cannibal(N::SpeciesInteractionNetwork)

    _cannibal = Any[]
    sp = species(N)

    for i in eachindex(sp)
        prey = collect(successors(N, sp[i]))

        if sp[i] ∈ prey
            push!(_cannibal, sp[i])
        end
    end

    return _cannibal
end

"""
Returns the percentage of species involved in a loop (motifs S3, D2, D4, D5, D6, D7)
"""
function loops(N::SpeciesInteractionNetwork)

    A = _get_matrix(N)

    # empty vector to puch spp index (proxy for id) to
    in_loop_spp = Any[]

    # look at loops that expand all the way to richness of network
    for i in 3:richness(N)
        # get the diagonal of the power transformed adj mat (indication of loops)
        _diag = diag(A^i)
        # get indices of spp with val > 0
        spp_ind = findall(!=(0), _diag)

        # add to master list
        if length(spp_ind) > 0
            append!(in_loop_spp, spp_ind)
        end
    end

    # return unique numbers (indices)
    return unique(in_loop_spp)
end

"""
remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Identifies and sets cannibalistic link to zero
"""
function remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    A = _get_matrix(N)

    A[diagind(A)] .= false

    nodes = Unipartite(species(N))
    edges = Binary(A)
    network = SpeciesInteractionNetwork(nodes, edges)

    return network
end

"""
max_sim(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Returns the mean maximum trophic similarity
"""
function max_sim(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    spp = species(N)
tsi = Any[]

for i in eachindex(spp)

    # we will update this with the actual max sim
    max_sim_i = 0.0

    # get preds and preys of spp i
    pred_i = successors(N, spp[i]) 
    prey_i = predecessors(N, spp[i])

    # cycle through the species in community and get their pred/prey to calc sim

    for j in eachindex(spp)
        # dont want to compare for same spp combo
        if j != i

            # get preds and preys of spp j
            pred_j = successors(N, spp[j]) 
            prey_j = predecessors(N, spp[j])

            shared_prey = length(collect(prey_i) ∩ collect(prey_j))
            shared_pred = length(collect(pred_i) ∩ collect(pred_j))
            all_prey = length(unique(vcat(collect(prey_i), collect(prey_j))))
            all_pred = length(unique(vcat(collect(pred_i), collect(pred_j))))

            similarity = (shared_prey + shared_pred)/(all_prey + all_pred)

            if similarity > max_sim_i
                max_sim_i = similarity
            end
        end
    end

    push!(tsi, max_sim_i)

end

    return mean(tsi)
end

"""
clustering(A::Matrix{Bool})

    Returns the mean clustering coefficient
"""
function clustering(A::Matrix{Bool})

    N = size(A, 1)
    
    # Calculate the Undirected Degree (k_i)
    # K_i is the total number of neighbors (in-degree + out-degree).
    A_undir = (A + A') .> 0 # A_undir[i,j] = 1 if there is a link i<->j or i->j or i<-j

    # The undirected degree k_i for species i is the sum of the i-th row (or column) of A_undir.
    k_undir = sum(A_undir, dims=2)[:]
    
    # Calculate the Number of Triangles (T_i)
    # In an undirected graph, the number of triangles T_i involving node i is half the (i, i) entry of A_undir^3.
    # We can calculate the total number of undirected links between neighbors of i directly.
    # The element (A_undir^2)_{ij} is the number of 2-paths between i and j.
    # The number of triangles T_i is the sum of links between the neighbors of i.
    
    # Let D be the number of cycles of length 3 (triangles)
    D = diag(A_undir^3) ./ 2

    # Calculate the Local Clustering Coefficient (C_i)
    C_values = Float64[] # Store local clustering coefficients

    for i in 1:N
        k_i = k_undir[i]
        
        # Denominator: Number of possible 2-paths (connections between neighbors)
        # This is the number of pairs of neighbors: k_i * (k_i - 1) / 2
        denominator = k_i * (k_i - 1) / 2
        
        if denominator == 0
            # Species with degree 0 or 1 cannot be part of a triangle.
            push!(C_values, 0.0) 
            continue
        end

        T_i = D[i] # Number of completed triangles involving node i
        
        # Local Clustering Coefficient C_i
        C_i = T_i / denominator
        push!(C_values, C_i)
    end
    
    # Calculate the Mean Clustering Coefficient
    mean_C = mean(C_values)
    
    return mean_C
end

"""
trophic_coherence(N::SpeciesInteractionNetwork)

Returns the trophic incoherence parameter q.
Lower q indicates higher trophic coherence.
"""
function trophic_coherence(N::SpeciesInteractionNetwork)

    A = _get_matrix(N)
    tl = trophic_level(N)

    spp = species(N)
    s = [tl[k] for k in spp]

    trophic_dist = Float64[]

    for i in eachindex(spp)
        for j in eachindex(spp)

            if A[i, j] == true
                push!(trophic_dist, s[i] - s[j])
            end

        end
    end

    # variance of trophic distances
    q = std(trophic_dist)

    return q
end

"""
trophic_variance(N::SpeciesInteractionNetwork)

Returns variance in trophic levels across species.
"""
function trophic_variance(N::SpeciesInteractionNetwork)

    tl = trophic_level(N)
    tls = collect(values(tl))

    return var(tls)

end

include("hopcroft_karp.jl")

"""
structural_controllability(N::SpeciesInteractionNetwork)

Returns the fraction of driver nodes required to control the network
following Liu et al. (2011). Lower values indicate higher controllability.
"""
function structural_controllability(N::SpeciesInteractionNetwork)
    A = _get_matrix(N)
    S = size(A, 1)

    # bipartite graph indices
    left  = collect(1:S)
    right = collect(S+1:2S)

    # Build adjacency list
    adj = Dict(u => Int[] for u in left)
    for i in 1:S
        for j in 1:S
            if A[i, j]
                push!(adj[i], S + j)
            end
        end
    end

    _, _, matching_size = hopcroft_karp_bipartite(left, right, adj)

    # Theoretical minimum: S - matching_size. 
    # If S == matching_size, you still need at least 1 driver node.
    Nd = max(1, S - matching_size)

    return Nd / S
end

function compute_reachable_to_top(N, top_set)

    reachable = Set(top_set)
    changed = true

    while changed
        changed = false

        for s in species(N)
            if any(n in reachable for n in predecessors(N, s))
                if !(s in reachable)
                    push!(reachable, s)
                    changed = true
                end
            end
        end
    end

    return reachable
end

# =========================
# CORE WORKFLOW
# =========================

function chain_metrics(N; max_depth=6)

    # --- STEP 1: Identify basal and top ---
    gen = SpeciesInteractionNetworks.generality(N)
    basal = collect(keys(filter(((k, v),) -> v == 0, gen)))

    vul = SpeciesInteractionNetworks.vulnerability(N)
    top_set = Set(keys(filter(((k, v),) -> v == 0, vul)))

    # If no structure exists → return early
    if isempty(basal) || isempty(top_set)
        return (ChLen = NaN, ChSD = NaN, ChNum = 0.0)
    end

    # --- STEP 2: Reachability pruning ---
    reachable = compute_reachable_to_top(N, top_set)

    # --- STEP 3: Memoized DFS ---
    memo = Dict{Any, Vector{Int}}()

    function dfs(node, visited, depth)

        if depth > max_depth
            return Int[]
        end

        if node in visited
            return Int[]
        end

        # prune unreachable nodes
        if node ∉ reachable
            return Int[]
        end

        if haskey(memo, node)
            return memo[node]
        end

        push!(visited, node)

        lengths = Int[]

        # If this is a top node → chain ends here
        if node in top_set
            push!(lengths, 0)
        end

        # Traverse UP the food web
        for nxt in predecessors(N, node)
            sub_lengths = dfs(nxt, visited, depth + 1)
            for l in sub_lengths
                push!(lengths, l + 1)
            end
        end

        delete!(visited, node)

        memo[node] = lengths
        return lengths
    end

    # --- STEP 4: Collect all chain lengths ---
    all_lengths = Int[]

    for b in basal
        append!(all_lengths, dfs(b, Set(), 0))
    end

    # --- STEP 5: Return summary stats ---
    if isempty(all_lengths)
        return (ChLen = NaN, ChSD = NaN, ChNum = 0.0)
    end

    return (
        ChLen = mean(all_lengths),
        ChSD  = std(all_lengths),
        ChNum = log(length(all_lengths))
    )
end
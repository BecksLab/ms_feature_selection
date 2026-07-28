# General sundry internal functions

using Extinctions
using FoodWebTools
using Graphs
using GraphsMatching
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics

# import other scripts with functions
include("hopcroft_karp.jl")

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
    tl = trophic_level(A)

    chain = chain_metrics(N; max_depth=6)

    # for centrality - freeman centralisation
    c = collect(values(centrality(N)))
    Cmax = maximum(c)

    # --- Robustness and Resilience with Error Handling ---
    n_reps = 100
    rob_vals = Float64[]
    res_vals = Float64[]
    
    attempts = 0
    while length(rob_vals) < n_reps && attempts < (n_reps * 2)
        attempts += 1
        try
            rb = robustness(N; threshold = 50, remove_disconnected = true)
            rs = resilience(extinction(N))
            
            push!(rob_vals, rb)
            push!(res_vals, rs)
        catch
            continue # Skip failed iterations
        end
    end

    D = Dict{Symbol,Any}(
        :richness => S,
        :links => L,
        :connectance => connectance(N),
        :diameter => FoodWebTools.diameter(A),
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
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)))/(richness(N)^2),
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)))/(richness(N)^2),
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)))/(richness(N)^2),
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)))/(richness(N)^2),
        :ρ => spectralradius(N),
        :centrality => sum(Cmax .- c),
        :loops => length(loops(N)) / S,
        :resilience => isempty(res_vals) ? NaN : mean(res_vals),
        :robustness => isempty(rob_vals) ? NaN : mean(rob_vals),
        :intervals => intervality(A),
        :MaxSim => max_sim(N),
        :Clust => clustering(A),
        :trophicCoherence => trophic_coherence(A),
        :trophicVar => trophic_variance(A),
        :control => structural_controllability(N),
        :Generality   => mean(gen),
        :Vulnerability => mean(vul),
        :MaxTL => maximum(values(tl))
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

_parser(x) = parse(Int, x)


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

    A = N.edges.edges
    sp = species(N)

    tl = trophic_level(A; species=sp)

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
trophic_variance(A::AbstractMatrix{Bool})

Returns variance in trophic levels across species.
"""
function trophic_variance(A::AbstractMatrix{Bool})

    tl = trophic_level(A)
    tls = collect(values(tl))

    return var(tls)

end

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
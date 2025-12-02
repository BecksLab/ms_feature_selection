# General sundry internal functions

using Extinctions
using Graphs
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics

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

    top2 = [k for (k, v) in _vul if v == 0]
    cl = [v for (k, v) in tl if k ∈ top2]
    if length(cl) == 0
        cl = 0.0
    end

    D = Dict{Symbol,Any}(
        :richness => S,
        :links => L,
        :connectance => connectance(N),
        :diameter => diameter(N),
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
        :ChLen => mean(cl),
        :ChSD => std(cl),
        :ChNum => log(length(cl)),
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
        :loops => length(loops(N)),
        :robustness => robustness(N; threshold = 50,
                                remove_disconnected = true))

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

    if richness(N) < 121

        # Create a directed graph (e.g., with richness of N)
        g = SimpleDiGraph(richness(N))

        # Build Di Graph
        # get matrix
        A = _get_matrix(N)

        # for every interaction pair that exists we push to DiGraph
        for x in axes(A, 1)
            for y in axes(A, 2)
                if A[x, y] == 1
                    add_edge!(g, x, y)
                end
            end   
        end

        # Find all simple cycles in the graph
        return simplecycles(g)
    else
        return Vector{Int64}[]
    end

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
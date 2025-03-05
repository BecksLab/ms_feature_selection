# General sundry internal functions

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
    gen = (collect(values(_gen)) * (1/(L/S)))
    vul = collect(values(SpeciesInteractionNetworks.vulnerability(N)))
    ind_maxgen = findmax(gen)[2]

    D = Dict{Symbol,Any}(
        :richness => S,
        :links => L,
        :connectance => connectance(N),
        :diameter => diameter(N),
        :complexity => complexity(N),
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :basal => sum(vec(sum(A, dims = 2) .== 0))/S,
        :top => sum(vec(sum(A, dims = 1) .== 0))/S,
        :generality => std(gen),
        :vulnerability => std(vul),
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], N))/S,
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], N))/S,
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], N))/S,
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], N))/S,
    )

    return D
end

"""
    maxrank(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}}

Returns the maximum possible rank of a Network
"""
function maxrank(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})
    return minimum(size(N))
end

"""
_get_matrix(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    Internal function to return a matrix of interactions from a
    SpeciesInteractionNetwork
"""
function _get_matrix(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    species = richness(N)
    n = zeros(Int64, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i, j] == true
                n[i, j] = 1
            end
        end
    end

    return n
end

"""
diameter(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Calculates the diameter of a food web. Where diameter is the longest 
    shortest path between two nodes
"""
function diameter(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # extract species names
    spp = species(N)
    # empty vector for storing shortest path for each spp
    shortpath = zeros(Int64, length(spp))

    # get shortest path
    for i in eachindex(spp)

        shortpath[i] = length(shortestpath(N, spp[i]))

    end

    # return max shortest path
    return findmax(shortpath)[1]
end

_parser(x) = parse(Int, x)

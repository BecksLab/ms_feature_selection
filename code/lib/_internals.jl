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
    gen = collect(values(_gen))
    _vul = SpeciesInteractionNetworks.vulnerability(N)
    vul = collect(values(_vul))
    ind_maxgen = findmax(gen)[2]
    l_s = L/S
    top = [k for (k,v) in _vul if v==0]
    tl = trophic_level(N)

    cl = [v for (k,v) in tl if k ∈ top]

    D = Dict{Symbol,Any}(
        :richness => S,
        :links => L,
        :connectance => connectance(N),
        :diameter => diameter(N),
        :complexity => complexity(N),
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :basal => sum(vec(sum(A, dims = 2) .== 0))/S,
        :top => sum(vec(sum(A, dims = 1) .== 0))/S,
        :l_S => l_s,
        :generality => std(gen)/l_s,
        :vulnerability => std(vul)/l_s,
        :trophic_level => mean(collect(values(tl))),
        :cl_mean => mean(cl),
        :cl_std => std(cl),
        :log_fc => log(length(cl)),
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], N)),
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], N)),
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], N)),
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], N)),
    )

    D[:intermediate] = 1 - D[:top] - D[:basal]

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

"""
trophic_level(N::SpeciesInteractionNetwork)

    Calculates the prey-averaged trophic level of all species in a network. 

    Williams, Richard J., and Neo D. Martinez. 2004. “Limits to Trophic Levels and Omnivory in Complex Food Webs: Theory and Data.” The American Naturalist 163 (3): 458–68. https://doi.org/10.1086/381964.
"""
function trophic_level(N::SpeciesInteractionNetwork)

    sp = species(N);


    # dictionary for path lengths
    pls = Dict{Any, Any}()

    for i in eachindex(sp) 
        # find shortest path to a basal species
        pls[sp[i]] = distancetobase(N, sp[i])
    end
    # return trophic level Dict
    return pls  
end
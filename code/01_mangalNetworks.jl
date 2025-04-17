using CSV
using DataFrames
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

import SpeciesInteractionNetworks.Mangal

include("lib/_internals.jl");

# query the ID number and other metadata for all ecological networks archived on mangal.io

number_of_networks = count(Mangal.MangalNetwork)
count_per_page = 100
number_of_pages = convert(Int, ceil(number_of_networks / count_per_page))

mangal_networks = DataFrame(
    id = Int64[],
    species = Int64[],
    links = Int64[],
    predators = Int64[],
    herbivores = Int64[],
    latitude = Any[],
    longitude = Any[],
    reference = Any[],
);

@showprogress "Paging networks" for page = 1:number_of_pages
    networks_in_page = Mangal.networks("count" => count_per_page, "page" => page - 1)
    @showprogress "Counting items" for current_network in networks_in_page
        D = Dict{Symbol,Any}()
        D[:id] = current_network.id
        D[:species] = count(Mangal.MangalNode, current_network)
        D[:links] = count(Mangal.MangalInteraction, current_network)
        D[:predators] =
            count(Mangal.MangalInteraction, current_network, "type" => "predation")
        D[:herbivores] =
            count(Mangal.MangalInteraction, current_network, "type" => "herbivory")
        if ismissing(current_network.position)
            D[:latitude] = string("NA")
            D[:longitude] = string("NA")
        else
            D[:latitude] = current_network.position[1]
            D[:longitude] = current_network.position[2]
        end
        D[:reference] = current_network.dataset.reference.bibtex
        push!(mangal_networks, D)
    end
end


## Filter for food webs (i.e. networks containing ONLY trophic links)
fw =
    (mangal_networks.predators .+ mangal_networks.herbivores) ./ mangal_networks.links .== 1
mangal_foodwebs = mangal_networks[fw, :]

# make a nice 'dataframe' to store network data
mangal_topology = DataFrame(
    id = Any[],
    richness = Any[],
    links = Any[],
    connectance = Any[],
    diameter = Any[],
    complexity = Any[],
    distance = Any[],
    basal = Any[],
    top = Any[],
    intermediate = Any[],
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
);

# make object to store each network so we can import it later
networks = DataFrame(id = Any[], network = Any[]);

@showprogress for i = 1:nrow(mangal_foodwebs)
    N = simplify(mangalnetwork(mangal_foodwebs.id[i]))

    N_d = Dict{Symbol,Any}()
    N_d[:id] = mangal_foodwebs.id[i]
    N_d[:network] = N
    push!(networks, N_d) # push network 'as is'

    N = simplify(render(Binary, N)) # make binary

    d = _network_summary(N)
    d[:id] = mangal_foodwebs.id[i]

    # send to df
    push!(mangal_topology, d)
end

# subset the initial networks query to only food webs used in `networks`
# this will act as the metadata
filter!(:id => x -> x ∈ networks.id, mangal_networks)

## Write networks as object
save_object("data/mangal/mangal_networks.jlds", networks)
## Write files
CSV.write("data/mangal/mangal_networks_metadata.csv", mangal_networks)
CSV.write("data/mangal/mangal_summary.csv", mangal_topology)

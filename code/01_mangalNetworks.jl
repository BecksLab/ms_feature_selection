using CSV
using DataFrames
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

import SpeciesInteractionNetworks.Mangal

# query the ID number and other metadata for all ecological networks archived on mangal.io

number_of_networks = count(Mangal.MangalNetwork)
count_per_page = 100
number_of_pages = convert(Int, ceil(number_of_networks/count_per_page))

mangal_networks = DataFrame(
    id = Int64[],
    species = Int64[],
    links = Int64[],
    predators = Int64[],
    herbivores = Int64[],
    latitude = Any[],
    longitude = Any[],
    );

@showprogress "Paging networks" for page in 1:number_of_pages
    networks_in_page = Mangal.networks("count" => count_per_page, "page" => page-1)
    @showprogress "Counting items" for current_network in networks_in_page
        D = Dict{Symbol,Any}()
        D[:id] = current_network.id
        D[:species] = count(Mangal.MangalNode, current_network)
        D[:links] = count(Mangal.MangalInteraction, current_network)
        D[:predators] = count(Mangal.MangalInteraction, current_network, "type" => "predation")
        D[:herbivores] = count(Mangal.MangalInteraction, current_network, "type" => "herbivory")
        if ismissing(current_network.position)
            D[:latitude] = string("NA")
            D[:longitude] = string("NA")
        else
            D[:latitude] = current_network.position[1]
            D[:longitude] = current_network.position[2]
        end
        push!(mangal_networks, D)
    end
end


## Filter for food webs (i.e. networks containing ONLY trophic links)
fw = (mangal_networks.predators .+ mangal_networks.herbivores) ./ mangal_networks.links .== 1
mangal_foodwebs = mangal_networks[fw, :]

# 20 datasets for now just for initial testing

# make a nice 'dataframe' to store network data
mangal_topology = DataFrame(
    id = Int64[],
    richness = Int64[],
    links = Int64[],
    connectance = Float64[],
    complexity = Float64[],
    distance = Float64[],
    basal = Float64[],
    top = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
);

# make object to store each network so we can import it later
networks = DataFrame(
    id = Any[],
    network = Any[]);

@showprogress for i in 1:nrow(mangal_foodwebs)
    N = simplify(mangalnetwork(mangal_foodwebs.id[i]))

        N_d = Dict{Symbol,Any}()
        N_d[:id] =mangal_foodwebs.id[i]
        N_d[:network] = N
        push!(networks, N_d) # push network 'as is'

        N = render(Binary, N) # make binary

        d = _network_summary(N)

        D = Dict{Symbol,Any}()
        D[:id] = mangal_foodwebs.id[i]
        D[:richness] = d[:richness]
        D[:links] = d[:links]
        D[:connectance] = d[:connectance]
        D[:complexity] = d[:complexity]
        D[:distance] = d[:distance]
        D[:basal] = d[:basal]
        D[:top] = d[:top]
        D[:S1] = d[:S1]
        D[:S2] = d[:S2]
        D[:S4] = d[:S4]
        D[:S5] = d[:S5]
        push!(mangal_topology, D)
end

## Write networks as object
save_object("data/raw/mangal/mangal_networks.jlds", networks)

## Write file
CSV.write("data/processed/mangal_networks_metadata.csv", mangal_foodwebs)
CSV.write("data/raw/mangal/mangal_networks.csv", networks)
CSV.write("data/processed/mangal_summary.csv", mangal_topology)
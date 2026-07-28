using CSV
using DataFrames
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

include("lib/_internals.jl");

# set seed
import Random
Random.seed!(66)

# get the name of all communities
web_names = readdir("data/vermaat_2009/raw")

# make empty dict to store network data
rows = Dict[];

# dataframe to store networks
networks = DataFrame(id = Any[], network = Any[]);

@showprogress "Getting network summary" for i in eachindex(web_names)

    web = web_names[i]

    # import data frame
    df = DataFrame(
        CSV.File(
            joinpath("data/vermaat_2009/raw", "$web"),
            header = ["consumer", "resource"],
            types = Symbol,
        ),
    )

    # get unique species (both consumer and resource)
    spp_list = unique(vcat(df.consumer, df.resource))

    nodes = Unipartite(spp_list)
    edges = Binary(zeros(Bool, length(spp_list), length(spp_list)))
    N = SpeciesInteractionNetwork(nodes, edges)

    # add interactions using interactions_all
    for j = 1:nrow(df)
        N[(df.consumer[j], df.resource[j])...] = true
    end

    # push network
    N_d = Dict{Symbol,Any}()
    N_d[:id] = replace.(web, ".csv" => "")
    N_d[:network] = N
    push!(networks, N_d)

    # push network summary
    d = _network_summary(N)
    d[:id] = replace.(web, ".csv" => "")

    # send to dict
    push!(rows, d)

end

# build DataFrame from rows dictionary
vermaat_topology = DataFrame(rows)

## Write networks as object
save_object("data/vermaat_2009/networks.jlds", networks)
## Write files
CSV.write("data/vermaat_2009/vermaat_summary.csv", vermaat_topology)

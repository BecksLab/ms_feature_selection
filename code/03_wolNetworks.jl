# weboflife_foodwebs.jl
using CSV
using DataFrames
using Downloads
using ZipFile
using SpeciesInteractionNetworks
using JLD2
using ProgressMeter
using JSON3
using StructTypes

include("lib/_internals.jl")

# Web of Life metadata endpoint
metadata_url = "https://www.web-of-life.es/get_networks.php"

println("Downloading Web of Life metadata...")

raw_json = Downloads.download(metadata_url)

json = JSON3.read(read(raw_json, String))

meta = DataFrame(json)

println("Loaded $(nrow(meta)) networks.")

# Filter food webs
foodweb_rows = filter(:network_name => x -> startswith(x, "FW_"), meta)

# unique networks
fw_ids = unique(foodweb_rows.network_name)

println("Found $(length(fw_ids)) food webs.")

# metadata dataframe

weboflife_networks = DataFrame(
    id = Any[],
    species = Any[],
    links = Any[],
    reference = Any[]
)

# prepare data frames
weboflife_topology = DataFrame(
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
    robustness = Any[],
    intervals = Any[],
    MaxSim = Any[],
    Clust = Any[],
    trophicCoherence = Any[],
    trophicVar = Any[],
    control = Any[]
)

networks = DataFrame(id = Any[], network = Any[])

# build networks

@showprogress "Building Web of Life networks" for fw in fw_ids

    df = filter(:network_name => x -> x == fw, foodweb_rows)

    predators = Symbol.(df.species2)
    prey = Symbol.(df.species1)

    spp = unique(vcat(prey, predators))

    A = zeros(Bool, length(spp), length(spp))

    sp_index = Dict(spp[i] => i for i in eachindex(spp))

    for r in eachrow(df)
        i = sp_index[Symbol(r.species1)]
        j = sp_index[Symbol(r.species2)]
        A[i,j] = 1
    end

    nodes = Unipartite(spp)
    edges = Binary(A)

    N = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)

    push!(networks, (id = fw, network = N))

    d = _network_summary(N)
    d[:id] = fw
    push!(weboflife_topology, d)

    push!(weboflife_networks, (
        id = fw,
        species = length(spp),
        links = sum(A),
        reference = "Web of Life database"
    ))
end

# Save results

save_object("data/weboflife/weboflife_networks.jlds", networks)
CSV.write(
    "data/weboflife/weboflife_networks_metadata.csv",
    weboflife_networks
)
CSV.write(
    "data/weboflife/weboflife_summary.csv",
    weboflife_topology
)
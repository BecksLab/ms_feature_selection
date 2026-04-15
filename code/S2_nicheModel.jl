using CSV
using DataFrames
using SpeciesInteractionNetworks
using ProgressMeter
using JLD2
import Random

include("lib/_internals.jl");

# catch time sinks
function safe_call(f; fallback = NaN)
    t = @async try f() catch; fallback end
    wait_time = 2.0  # seconds

    for _ in 1:Int(wait_time * 10)
        if istaskdone(t)
            return fetch(t)
        end
        sleep(0.1)
    end

    return fallback  # timeout
end

Random.seed!(66)

net_summs = CSV.read("data/cleaned/all_networks.csv", DataFrame)

n_reps = 100

# collect rows first
rows = Vector{Dict{Symbol,Any}}()

total = nrow(net_summs) * n_reps
p = Progress(total)

for j in 1:nrow(net_summs)

    S_input = net_summs.richness[j]
    C = net_summs.connectance[j]

    for i in 1:n_reps
        next!(p)

        try
            # Build network
            N = structuralmodel(NicheModel, S_input, C)

            L = links(N)
            S = richness(N)

            # basic guards
            if L == 0 || S == 0
                continue
            end

            l_s = L / S
            if l_s == 0
                continue
            end

            # adjacency
            A = _get_matrix(N)

            # generality & vulnerability
            _gen = SpeciesInteractionNetworks.generality(N)
            _vul = SpeciesInteractionNetworks.vulnerability(N)

            gen = collect(values(_gen))
            vul = collect(values(_vul))

            if isempty(gen) || isempty(vul)
                continue
            end

            # safe indexing
            ind_maxgen = try findmax(gen)[2] catch; continue end
            species_keys = collect(keys(_gen))

            if ind_maxgen > length(species_keys)
                continue
            end

            # trophic structure
            top = sum(vec(sum(A, dims = 1) .== 0))
            basal = sum(vec(sum(A, dims = 2) .== 0))
            int = S - (basal + top)

            # safe metrics
            diam = try diameter(A) catch; NaN end
            dist = try distancetobase(N, species_keys[ind_maxgen]) catch; NaN end
            tl = try trophic_level(N) catch; Dict() end
            tl_mean = isempty(tl) ? NaN : mean(Float64.(values(tl)))

            chain = safe_call(() -> chain_metrics(N; max_depth = 4), fallback = nothing)

            path_mean = try mean(Float64.(pathlengths(N))) catch; NaN end

            deg = try values(SpeciesInteractionNetworks.degree(N)) catch; [] end
            linksd = isempty(deg) ? NaN : std(deg) / l_s

            predprey = (top + int) == 0 ? NaN : (basal + int)/(top + int)

            # for centrality - freeman centralisation
            c = collect(values(centrality(N)))
            Cmax = maximum(c)

            d2 = Dict{Symbol,Any}(
                :id => j,
                :run => i,
                :richness => S,
                :links => L,
                :connectance => connectance(N),
                :diameter => diam,
                :distance => dist,
                :basal => basal / S,
                :top => top / S,
                :intermediate => int / S,
                :herbivory => length(herbivore(N)) / S,
                :omnivory => length(omnivore(N)) / S,
                :cannibal => length(cannibal(N)) / S,
                :predpreyRatio => predprey,
                :l_S => l_s,
                :GenSD => std(gen) / l_s,
                :VulSD => std(vul) / l_s,
                :TL => tl_mean,
                :ChLen => chain === nothing ? NaN : chain.ChLen,
                :ChSD  => chain === nothing ? NaN : chain.ChSD,
                :ChNum => chain === nothing ? NaN : chain.ChNum,
                :path => path_mean,
                :LinkSD => linksd,
                :S1 => try length(findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)))/(S^2) catch; NaN end,
                :S2 => try length(findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)))/(S^2) catch; NaN end,
                :S4 => try length(findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)))/(S^2) catch; NaN end,
                :S5 => try length(findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)))/(S^2) catch; NaN end,
                :centrality => try sum(Cmax .- c) catch; NaN end,
                :loops => try length(loops(N)) / S catch; NaN end,
                :intervals => try intervality(A) catch; NaN end,
                :MaxSim => try max_sim(N) catch; NaN end,
                :Clust => try clustering(A) catch; NaN end,
                :trophicCoherence => try trophic_coherence(N) catch; NaN end,
                :trophicVar => try trophic_variance(N) catch; NaN end,
            )

            push!(rows, d2)

        catch e
            println("Error at j=$j, i=$i")
            println(e)
            continue
        end
    end
end

# build DataFrame
niche_topology = DataFrame(rows)

CSV.write("data/cleaned/nicheNetworks.csv", niche_topology)
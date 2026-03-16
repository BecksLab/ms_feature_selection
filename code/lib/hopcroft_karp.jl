"""
    hopcroft_karp_bipartite(left, right, adj)

Hopcroft–Karp on bipartite graph defined by:

- `left`:  vector of left node ids
- `right`: vector of right node ids
- `adj`:   adjacency list: Dict{Int,Vector{Int}} mapping left -> neighbors in right

Returns:
- `pairU`: Dict mapping each left to its matched right (or 0)
- `pairV`: Dict mapping each right to its matched left (or 0)
- `matching_size`: Int
"""
function hopcroft_karp_bipartite(left, right, adj)
    # All nodes include 0 as NIL
    INF = typemax(Int)

    # match maps
    pairU = Dict(u => 0 for u in left)
    pairV = Dict(v => 0 for v in right)

    dist = Dict(u => 0 for u in left)

    bfs() = begin
        queue = []
        for u in left
            if pairU[u] == 0
                dist[u] = 0
                push!(queue, u)
            else
                dist[u] = INF
            end
        end
        dist[0] = INF

        while !isempty(queue)
            u = popfirst!(queue)
            if dist[u] < dist[0]
                for v in adj[u]
                    if dist[pairV[v]] == INF
                        dist[pairV[v]] = dist[u] + 1
                        push!(queue, pairV[v])
                    end
                end
            end
        end
        return dist[0] != INF
    end

    dfs(u) = begin
        if u != 0
            for v in adj[u]
                if dist[pairV[v]] == dist[u] + 1 &&
                   dfs(pairV[v])
                    pairV[v] = u
                    pairU[u] = v
                    return true
                end
            end
            dist[u] = INF
            return false
        end
        return true
    end

    matching = 0
    while bfs()
        for u in left
            if pairU[u] == 0 && dfs(u)
                matching += 1
            end
        end
    end

    return pairU, pairV, matching
end
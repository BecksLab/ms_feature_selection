using Graphs

"""
    diameter(A::AbstractMatrix{<:Bool})

Calculates the diameter of a food web represented by an adjacency matrix A.
It explicitly builds the graph using the edge list (A[i, j] = 1 means i -> j)
to ensure correct construction, then uses the Floyd-Warshall algorithm from Graphs.jl.
The diameter is filtered using the V-1 rule to exclude unreachable paths.

"""
function diameter(A::AbstractMatrix{<:Bool})
    # 1. Explicitly build the graph using an edge list (most robust construction)
    num_vertices = size(A, 1)
    
    food_web_graph = SimpleDiGraph(num_vertices)
    
    # Add edges based on the standard convention: A[i, j] = 1 means i -> j
    for i in 1:num_vertices
        for j in 1:num_vertices
            if A[i, j] == 1
                add_edge!(food_web_graph, i, j)
            end
        end
    end
    
    # 2. Calculate All-Pairs Shortest Paths
    dists = floyd_warshall_shortest_paths(food_web_graph)
    distance_matrix = dists.dists

    # 3. Calculate the Diameter using the V-1 filter
    # The longest possible finite shortest path is V - 1.
    max_finite_path_length = num_vertices - 1 

    food_web_diameter = 0
    
    for d in distance_matrix
        # Filters out self-loops (d=0) and unreachable paths (d > max_finite_path_length)
        if d > 0 && d <= max_finite_path_length
            if d > food_web_diameter
                food_web_diameter = d
            end
        end
    end

    return food_web_diameter
end

# Example 1: Long Trophic Chain (Diameter = 3)
A1 = [
    0 1 0 0  # 1 eats 2
    0 0 1 0  # 2 eats 3
    0 0 0 1  # 3 eats 4
    0 0 0 0
]

diameter1 = diameter(A1)

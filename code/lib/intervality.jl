using LinearAlgebra # Required for size()

# --- Helper Function for Efficient State Swapping ---
"""
swap_matrix_rows_cols!(A, i, j)

Swaps row i with row j, and column i with column j, in place.
This is necessary to maintain the adjacency matrix relationship after a permutation.
"""
function swap_matrix_rows_cols!(A::Matrix{Bool}, i::Int, j::Int)
    # Check bounds
    if i == j || i < 1 || j < 1 || i > size(A, 1) || j > size(A, 1)
        return A
    end
    
    # Swap Rows i and j
    A[i, :], A[j, :] = A[j, :], A[i, :]
    # Swap Columns i and j
    A[:, i], A[:, j] = A[:, j], A[:, i]
    return A
end

# --- Core Function: Calculate Gaps ---
"""
calculate_gaps(A::Matrix{Int})

Calculates the total number of diet gaps (G) for a given ordered adjacency matrix A.
Gaps occur when a consumer eats two resources but skips an intermediate resource 
in the current niche axis ordering.
"""
function calculate_gaps(A::Matrix{Bool})
    S = size(A, 1)
    total_gaps = 0

    # Iterate over consumers (rows)
    for i in 1:S
        # Find the column indices (positions) of the resources consumed
        resources_consumed = findall(A[i, :] .== 1)

        if !isempty(resources_consumed)
            # k_min is the position of the lowest-ranked prey
            k_min = minimum(resources_consumed)
            # k_max is the position of the highest-ranked prey
            k_max = maximum(resources_consumed)

            # Potential diet size: the length of the continuous block spanned by the diet
            potential_diet_size = k_max - k_min + 1

            # Actual diet size: the number of ones in the row
            actual_diet_size = length(resources_consumed)

            # Gaps = (Potential Size - Actual Size)
            total_gaps += (potential_diet_size - actual_diet_size)
        end
    end
    return total_gaps
end

# --- Main Function: Simulated Annealing Optimization ---
"""
intervality(A_initial::Matrix{Bool}; 
                         T_initial=10.0, cooling_rate=0.99999, iterations=500000)

Uses Simulated Annealing to find the minimum number of gaps (G_min) by searching
for the optimal species ordering (permutation) of the adjacency matrix.
"""
function intervality(A_initial::Matrix{Bool}; 
                                  T_initial=10.0, 
                                  cooling_rate=0.99999, 
                                  iterations=500000)
    
    S = size(A_initial, 1)
    
    # 1. Initial State Setup
    current_A = copy(A_initial)
    current_gaps = calculate_gaps(current_A)
    best_gaps = current_gaps
    
    T = T_initial 

    # 2. SA Optimization Loop
    for k in 1:iterations
        # 3. Propose New State: Swap two random positions (indices)
        idx1, idx2 = rand(1:S, 2)
        while idx1 == idx2
            idx2 = rand(1:S)
        end
        
        # Create new state by copying and swapping the matrix
        new_A = copy(current_A)
        swap_matrix_rows_cols!(new_A, idx1, idx2)

        new_gaps = calculate_gaps(new_A)
        
        # 4. Acceptance Check (Metropolis Criterion)
        delta_gaps = new_gaps - current_gaps
        
        if delta_gaps <= 0  # Accept if better or equal
            current_gaps = new_gaps
            current_A = new_A 
            
            # Update the overall best solution
            if new_gaps < best_gaps
                best_gaps = new_gaps
            end
            
        else # Accept a worse state with probability p
            acceptance_prob = exp(-delta_gaps / T)
            if rand() < acceptance_prob
                current_gaps = new_gaps
                current_A = new_A 
            end
        end
        
        # 5. Cooling Schedule (Geometric)
        T *= cooling_rate
        
        # Early exit (if perfect interval graph is found)
        if best_gaps == 0
            break
        end
    end

    return best_gaps
end
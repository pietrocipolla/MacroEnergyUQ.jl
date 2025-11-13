"""
    transport_ot(source_marginal::Vector{Float64}, target_marginal::Vector{Float64}, cost_matrix::Matrix{Float64})

Compute optimal transport plan using the R package 'transport'.

This function solves the optimal transport problem between two probability distributions
using the network simplex algorithm from the R transport package.

# Arguments
- `source_marginal`: Vector of probability weights for the source distribution (must sum to 1)
- `target_marginal`: Vector of probability weights for the target distribution (must sum to 1)
- `cost_matrix`: Matrix of pairwise costs between source and target points (size: m x n)

# Returns
- Vector of integers where result[i] indicates which target point source point i is mapped to

# Example
```julia
n_data = 100
n_centroids = 4
data_marginal = fill(1/n_data, n_data)
centroids_marginal = fill(1/n_centroids, n_centroids)
C = pairwise(SqEuclidean(), data, centroids; dims=2)
cluster_assignments = transport_ot(data_marginal, centroids_marginal, C)
```
"""
function transport_ot(source_marginal::Vector{Float64}, target_marginal::Vector{Float64}, cost_matrix::Matrix{Float64})
    m = length(source_marginal)
    n = length(target_marginal)
    
    # Verify dimensions
    if size(cost_matrix) != (m, n)
        throw(ArgumentError("Cost matrix must be m x n where m = length(source_marginal) and n = length(target_marginal)"))
    end
    
    if !isapprox(sum(source_marginal), 1.0, atol=1e-10)
        throw(ArgumentError("Source marginal must sum to 1"))
    end
    
    if !isapprox(sum(target_marginal), 1.0, atol=1e-10)
        throw(ArgumentError("Target marginal must sum to 1"))
    end
    
    # Transfer data to R
    @rput source_marginal
    @rput target_marginal
    @rput cost_matrix
    
    # Solve optimal transport using R's transport package
    R"""
    # Load transport package (install if not available)
    if (!require("transport", quietly = TRUE)) {
        install.packages("transport", repos = "https://cloud.r-project.org")
        library(transport)
    }
    
    # Solve optimal transport
    transport_plan <- transport(source_marginal, target_marginal, cost_matrix)
    
    # Create assignment vector: for each source point, assign it to the target with maximum mass
    # Initialize with zeros
    assignment <- integer(length(source_marginal))
    
    # For each source point, find the target with maximum transported mass
    for (i in unique(transport_plan$from)) {
        # Get all targets for this source
        targets <- transport_plan[transport_plan$from == i, ]
        # Assign to target with maximum mass
        assignment[i] <- targets$to[which.max(targets$mass)]
    }
    """
    
    # Get the assignment back from R
    assignment = rcopy(R"as.integer(assignment)")
    
    return assignment
end

"""
    transport_ranking(marginal::Vector{Float64}, cost_matrix::Matrix{Float64})

Compute optimal transport ranking using the R package 'transport'.

This is a convenience wrapper for transport_ot when source and target distributions are identical.

# Arguments
- `marginal`: Vector of probability weights for both source and target distributions (must sum to 1)
- `cost_matrix`: Matrix of pairwise costs between source and target points (size: n x n)

# Returns
- Vector of integers where result[i] indicates which target point source point i is mapped to

# Example
```julia
n = 100
marginal = fill(1/n, n)
cost_matrix = pairwise(SqEuclidean(), data, grid; dims=2)
ranking = transport_ranking(marginal, cost_matrix)
```
"""
function transport_ranking(marginal::Vector{Float64}, cost_matrix::Matrix{Float64})
    return transport_ot(marginal, marginal, cost_matrix)
end

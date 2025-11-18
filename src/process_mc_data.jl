"""
    process_mc_data(data::Matrix{Float64}, n_threads::Int, starting_point::Vector{Float64})

Prepare a matrix of input samples for model run. It splits the data into `n_threads` chunks and reorder each chunk from the closest to the farthest from the center of the parameter space defined by `starting_point`.

Arguments:
- `data`: Data matrix.
- `n_threads`: Number of parallel threads to use
- `quantile_transform`: Boolean indicating whether to transform the data to [0,1]^d space using Wasserstein rankings. The transformation preserves the multivariate structure of the data while ensuring uniform marginals. Computationally expensive for large datasets.
- `sorting_algorithm`: Symbol indicating the sorting algorithm to use (currently :nearest_neighbor and :tsp are supported)
- `starting_point`: Parameter vector that defines the starting point for the nearest_neighbor reordering algorithm (must have the same length as the number of columns in data)
- `optimizer`: An optimizer to use when `sorting_algorithm=:tsp`

Throws ArgumentError if:
- number of columns in data doesn't match the length of starting_point
- n_threads is less than 1

Returns:
- A matrix with same values of `data` but with reordered rows, and an additional column that represents the cluster for parallelization.
"""
function process_mc_data(data::Matrix{Float64}, 
                        n_threads::Int; 
                        quantile_transform::Bool = true,
                        sorting_algorithm::Symbol = :nearest_neighbor,
                        starting_point::Vector{Float64} = fill(0.5, size(data, 1)), 
                        optimizer=nothing)
    # Check the sizes
    if size(data, 1) != length(starting_point)
        throw(ArgumentError("Number of columns in the matrix ($(size(data, 1))) must match the length of the parameter vector ($(length(starting_point)))"))
    end
    
    if n_threads < 1
        throw(ArgumentError("Number of threads must be at least 1"))
    end

    if n_threads > size(data, 2)
        throw(ArgumentError("Number of threads ($n_threads) is greater than number of points ($(size(data, 2)))."))
    end

    if sorting_algorithm != :nearest_neighbor && sorting_algorithm != :tsp
        throw(ArgumentError("Unsupported sorting algorithm: $sorting_algorithm. Currently, only :nearest_neighbor and :tsp are supported."))
    end

    if sorting_algorithm == :tsp && isnothing(optimizer)
        throw(ArgumentError("An optimizer must be provided when using the :tsp sorting algorithm."))
    end
    
    # Number of columns (data points) in the data matrix
    n_cols = size(data, 2)
    
    # Transform data to [0,1]^d space using Wasserstein rankings
    data_quantiles, starting_point_quantiles = data, starting_point
    if quantile_transform
        data_quantiles, starting_point_quantiles = _multivariate_quantiles(data, starting_point)
    end

    # Generate the clusters for parallelization
    if n_threads == 1
        clusters = ones(Int, n_cols)
    else
        clusters = _generate_cluster(data_quantiles, n_threads)
    end

    # Reorder points within each cluster
    data, original_indices = _reorder_within_clusters(data, 
                                                      data_quantiles, 
                                                      clusters, 
                                                      starting_point_quantiles, 
                                                      sorting_algorithm,
                                                      optimizer)

    return data, clusters, original_indices
end

"""
    _multivariate_quantiles(data::Matrix{Float64}, starting_point::Vector{Float64})

Private function that transforms data points to uniform quantiles in [0,1]^d space using Optimal Transport.
The transformation preserves the multivariate structure of the data while ensuring uniform marginals.

The function:
1. Combines the input data and starting_point into a single matrix
2. Generates a uniform grid using Sobol sequences as the target distribution
3. Computes the optimal transport plan between the empirical data distribution and the uniform grid
4. Uses the transport plan to map each point (including starting_point) to its corresponding quantile position

Arguments:
- `data`: Matrix where each column represents a point in the original space
- `starting_point`: Vector representing a reference point in the original space

Returns:
- Tuple containing:
  1. Matrix of transformed data points in [0,1]^d space
  2. Vector of transformed starting_point in [0,1]^d space

Note: The transformation is based on the Earth Mover's Distance (EMD) optimal transport solution
"""
function _multivariate_quantiles(data::Matrix{Float64}, starting_point::Vector{Float64})
    # Join data and starting_point
    data = hcat(data, starting_point)

    n = size(data, 2)
    d = size(data, 1)
    
    # Generate a uniform grid using Sobol sequences
    grid = QuasiMonteCarlo.sample(n, d, SobolSample())
    
    # Compute the cost matrix
    C = pairwise(SqEuclidean(), data, grid; dims=2)
    
    # Define the marginal measures for the OT problem
    marginal = fill(1 / n, n)
    
    # Solve the OT problem to get the transport plan
    ranking = transport_ranking(marginal, C)

    # Compute the quantiles for each dimension
    quantile_points = grid[:, ranking[1:(end-1)]]
    starting_point_quantiles = grid[:, ranking[end]]
    
    return quantile_points, starting_point_quantiles    
end

"""
    _generate_cluster(data::Matrix{Float64}, n_threads::Int)

Private function that generates cluster assignments using Optimal Transport (OT).

The function performs the following steps:
1. Normalizes the data to [0,1]^d space using Wasserstein rankings:
   - Generates a uniform grid using Sobol sequences
   - Computes optimal transport between data and grid points
   - Uses the transport plan to map data points to normalized positions
2. Generates `n_threads` centroids using Sobol sequences
3. Solves another OT problem to assign points to clusters:
   - Source: normalized data points with uniform distribution
   - Target: centroids with uniform distribution
   - Cost: squared Euclidean distances
4. Assigns each point to the centroid assigned by the optimal transport

Arguments:
- `data`: Input data matrix where each column is a point in the parameter space
- `n_threads`: Number of clusters to generate (matches desired number of parallel threads)

Returns:
- A row vector of cluster assignments where each element is an integer in [1, n_threads]
"""
function _generate_cluster(data::Matrix{Float64}, n_threads::Int)
    n_rows = size(data, 1)
    n_cols = size(data, 2)

    # Cluster the normalized data using OT
    # Define the centroids
    centroids = QuasiMonteCarlo.sample(n_threads, n_rows, SobolSample())

    # Define the marginal measures for the OT problem
    data_marginal = fill(1 / n_cols, n_cols)
    centroids_marginal = fill(1 / n_threads, n_threads)

    # Compute the cost matrix
    C = pairwise(SqEuclidean(), data, centroids; dims=2)

    # Solve OT problem using R transport package
    cluster_assignments = transport_ot(data_marginal, centroids_marginal, C)

    return cluster_assignments
end

"""
    _reorder_within_clusters(data::Matrix{Float64}, data_quantiles::Matrix{Float64}, clusters::Vector{Int}, starting_point_quantile::Vector{Float64}, sorting_algorithm::Symbol)

Private function that reorders points within each cluster.
If `sorting_algorithm=:nearest_neighbor`, it builds a path starting from the point closest to starting_point_quantile by always moving to the nearest unvisited point.
If `sorting_algorithm=:tsp`, it uses a traveling salesperson problem approach to find the optimal path within each cluster.

Arguments:
- `data`: Matrix of original input data
- `data_quantiles`: Matrix with the multivariate quantiles
- `clusters`: Vector of cluster assignments for each point
- `starting_point_quantile`: Reference point in quantile space for starting the path
- `sorting_algorithm`: Symbol indicating the sorting algorithm to use (currently :nearest_neighbor and :tsp are supported)

Returns:
- Tuple containing:
  1. Reordered original data matrix
  2. Indices giving the new order of points
"""
function _reorder_within_clusters(data::Matrix{Float64}, 
                                data_quantiles::Matrix{Float64}, 
                                clusters::Vector{Int},
                                starting_point_quantile::Vector{Float64},
                                sorting_algorithm::Symbol,
                                optimizer)
    
    if sorting_algorithm == :nearest_neighbor
        return _reorder_within_clusters_nn(data, data_quantiles, clusters, starting_point_quantile)
    elseif sorting_algorithm == :tsp
        return _reorder_within_clusters_tsp(data, data_quantiles, clusters, optimizer)
    end
end

"""
    _reorder_within_clusters_nn(data::Matrix{Float64}, data_quantiles::Matrix{Float64}, clusters::Vector{Int}, starting_point_quantile::Vector{Float64})

Private function that reorders points within each cluster using a nearest-neighbor approach.
Starting from the point closest to starting_point_quantile, it builds a path by always moving to the nearest unvisited point.

Arguments:
- `data`: Matrix of original input data
- `data_quantiles`: Matrix with the multivariate quantiles
- `clusters`: Vector of cluster assignments for each point
- `starting_point_quantile`: Reference point in quantile space for starting the path

Returns:
- Tuple containing:
  1. Reordered original data matrix
  2. Indices giving the new order of points
"""
function _reorder_within_clusters_nn(data::Matrix{Float64}, 
                                     data_quantiles::Matrix{Float64}, 
                                     clusters::Vector{Int},
                                     starting_point_quantile::Vector{Float64})
    n_clusters = maximum(clusters)
    n_cols = size(data, 2)
    
    # Initialize output arrays
    reordered_data = zeros(size(data))
    reordered_indices = collect(1:n_cols)
    
    # Process each cluster
    for cluster in 1:n_clusters
        # Find points belonging to current cluster
        cluster_mask = clusters .== cluster
        cluster_points = data_quantiles[:, cluster_mask]
        original_indices = findall(cluster_mask)
        n_points = length(original_indices)
        
        if n_points > 0
            # Array to store the new order of points
            new_order = Vector{Int}(undef, n_points)
            
            # Find the point closest to starting_point_quantile
            distances_to_starting_point = [sum((cluster_points[:, i] .- starting_point_quantile).^2) for i in 1:n_points]
            current_idx = argmin(distances_to_starting_point)
            new_order[1] = current_idx
            
            # Build path using nearest neighbor
            unvisited = trues(n_points)
            unvisited[current_idx] = false
            
            # Find remaining points using nearest neighbor
            for i in 2:n_points
                current_point = cluster_points[:, current_idx]
                
                # Calculate distances to all unvisited points
                distances = [unvisited[j] ? sum((cluster_points[:, j] .- current_point).^2) : Inf 
                           for j in 1:n_points]
                
                # Find nearest unvisited point
                current_idx = argmin(distances)
                new_order[i] = current_idx
                unvisited[current_idx] = false
            end
            
            # Apply the new ordering to both data matrices
            reordered_data[:, original_indices] = data[:, original_indices[new_order]]
            reordered_indices[original_indices] = original_indices[new_order]
        end
    end
    
    return reordered_data, reordered_indices
end

"""
    _reorder_within_clusters_tsp(data::Matrix{Float64}, data_quantiles::Matrix{Float64}, clusters::Vector{Int}, starting_point_quantile::Vector{Float64})

Private function that reorders points within each cluster using a traveling salesperson problem approach.
Starting from the point closest to starting_point_quantile, it builds a path by always moving to the nearest unvisited point.

Arguments:
- `data`: Matrix of original input data
- `data_quantiles`: Matrix with the multivariate quantiles
- `clusters`: Vector of cluster assignments for each point
- `starting_point_quantile`: Reference point in quantile space for starting the path

Returns:
- Tuple containing:
  1. Reordered original data matrix
  2. Indices giving the new order of points
"""
function _reorder_within_clusters_tsp(data::Matrix{Float64}, 
                                      data_quantiles::Matrix{Float64}, 
                                      clusters::Vector{Int},
                                      optimizer)
    n_clusters = maximum(clusters)
    n_cols = size(data, 2)

    # Compute the cost matrix
    C = pairwise(SqEuclidean(), data_quantiles, data_quantiles; dims=2);
    # Add a zero row and column at the beginning. This should allow to avoid closed-loop tours
    # Suggested in https://stackoverflow.com/questions/6733999/what-is-the-problem-name-for-traveling-salesman-problemtsp-without-considering
    C = [0 zeros(1, size(C, 2));
        zeros(size(C, 1), 1) C]
    
    # Initialize output arrays
    reordered_data = zeros(size(data))
    reordered_indices = collect(1:n_cols)
    
    # Process each cluster
    for cluster in 1:n_clusters
        # Find points belonging to current cluster
        cluster_mask = vcat(true, clusters .== cluster)
        cluster_C = C[cluster_mask, cluster_mask]
        original_indices = findall(cluster_mask)[2:end] .- 1  # Adjust for added zero row/column
        n_points = length(original_indices) + 1
        
        if n_points > 0
            iterative_model = build_tsp_model(cluster_C, n_points, optimizer)
            optimize!(iterative_model)
            assert_is_solved_and_feasible(iterative_model)
            cycle = subtour(iterative_model[:x])
            
            while 1 < length(cycle) < n_points
                S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
                @constraint(
                    iterative_model,
                    sum(iterative_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
                )
                optimize!(iterative_model)
                assert_is_solved_and_feasible(iterative_model)
                cycle = subtour(iterative_model[:x])
            end

            # Extract the tour order adjusting for added zero row/column
            new_order = get_tour_order(iterative_model) .- 1 

            # Apply the new ordering to both data matrices
            reordered_data[:, original_indices] = data[:, original_indices[new_order]]
            reordered_indices[original_indices] = original_indices[new_order]
        end
    end
    
    return reordered_data, reordered_indices
end


# TSP functions from https://jump.dev/JuMP.jl/stable/tutorials/algorithms/tsp_lazy_constraints/
function build_tsp_model(d, n, optimizer)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @objective(model, Min, sum(d .* x) / 2)
    @constraint(model, [i in 1:n], sum(x[i, :]) == 2)
    @constraint(model, [i in 1:n], x[i, i] == 0)
    return model
end

function subtour(edges::Vector{Tuple{Int,Int}}, n)
    shortest_subtour, unvisited = collect(1:n), Set(collect(1:n))
    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors =
                [j for (i, j) in edges if i == current && j in unvisited]
        end
        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
        end
    end
    return shortest_subtour
end

function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

"""
    get_tour_order(model::JuMP.Model)

Extract the tour order from a solved JuMP TSP model.

# Arguments
- `model`: A solved JuMP model.

# Returns
A vector of node indices representing the tour order.
"""
function get_tour_order(model::Model)
    # If x is an array (as in the JuMP TSP example)
    n = size(model[:x], 1)
    xval = value.(model[:x])
    
    # Build adjacency dictionary
    connections = Dict(i => findall(j -> xval[i, j] > 0.5, 1:n) for i in 1:n)

    # Reconstruct the tour order
    current = connections[1][1]
    tour = [1, current]  # start at point 1 (arbitrary)
    precedent = 1
    
    for _ in 1:n-1
        candidates = connections[current]
        new_current = candidates[candidates .!= precedent][1]
        precedent = current
        current = new_current
        push!(tour, current)
    end
    
    return tour[2:(end-1)]  # Exclude the added starting point
end
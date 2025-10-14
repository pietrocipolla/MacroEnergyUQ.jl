"""
    process_mc_data(data::Matrix{Float64}, n_threads::Int, params::Vector{Float64})

Prepare a matrix of input samples for model run. It splits the data into `n_threads` chunks and reorder each chunk from the closest to the farthest from the center of the parameter space defined by `params`.

Arguments:
- `data`: Data matrix.
- `n_threads`: Number of parallel threads to use
- `params`: Parameter vector that defines the starting point for the reordering algorithm (must have the same length as the number of columns in data)

Throws ArgumentError if:
- number of columns in data doesn't match the length of params
- n_threads is less than 1

Returns:
- A matrix with same values of `data` but with reordered rows, and an additional column that represents the cluster for parallelization.
"""
function process_mc_data(data::Matrix{Float64}, n_threads::Int, params::Vector{Float64})
    # Check the sizes
    if size(data, 1) != length(params)
        throw(ArgumentError("Number of columns in the matrix ($(size(data, 1))) must match the length of the parameter vector ($(length(params)))"))
    end
    
    if n_threads < 1
        throw(ArgumentError("Number of threads must be at least 1"))
    end
    
    # Number of rows in the data matrix
    n_rows = size(data, 1)
    n_cols = size(data, 2)
    
    # Transform data to [0,1]^d space using Wasserstein rankings
    data_quantiles, params_quantiles = _multivariate_quantiles(data, params)

    # Generate the clusters for parallelization
    if n_threads == 1
        cluster_assignments = ones(Int, n_cols)
    else
        cluster_assignments = _generate_cluster(data_quantiles, n_threads)
    end

    # Reorder points within each cluster
    data = _reorder_within_clusters(data, data_quantiles, params_quantiles)

    return data, cluster_assignments
end

function _multivariate_quantiles(data::Matrix{Float64}, params::Vector{Float64})
    # Join data and params
    data = hcat(data, params)

    n = size(data, 2)
    d = size(data, 1)
    
    # Generate a uniform grid using Sobol sequences
    grid = QuasiMonteCarlo.sample(n, d, SobolSample())
    
    # Compute the cost matrix
    C = pairwise(SqEuclidean(), data, grid; dims=2)
    
    # Define the marginal measures for the OT problem
    marginal = fill(1 / n, n)
    
    # Solve the OT problem to get the transport plan
    transport_plan = PythonOT.emd(marginal, marginal, C)
    ranking = [argmax(transport_plan[i, :]) for i in 1:n]
    
    # Compute the quantiles for each dimension
    quantile_points = grid[:, ranking[1:(end-1)]]
    params_quantiles = grid[:, ranking[end]]
    
    return quantile_points, params_quantiles    
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
4. Assigns each point to the centroid with maximum transport probability

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
    C = pairwise(SqEuclidean(), data, centroids; dims=2);

    ot_plan = PythonOT.emd(data_marginal, centroids_marginal, C);

    # Find the cluster assignment for each point (index of maximum value in each row)
    cluster_assignments = [argmax(ot_plan[i,:]) for i in 1:n_cols]

    return cluster_assignments'
end

"""
    _reorder_within_clusters(data::Matrix{Float64}, data_quantiles::Matrix{Float64}, clusters::Vector{Int}, params_quantile::Vector{Float64})

Private function that reorders points within each cluster using a nearest-neighbor approach.
Starting from the point closest to params_quantile, it builds a path by always moving to the nearest unvisited point.

Arguments:
- `data`: Matrix of original input data
- `data_quantiles`: Matrix with the multivariate quantiles
- `clusters`: Vector of cluster assignments for each point
- `params_quantile`: Reference point in quantile space for starting the path

Returns:
- Tuple containing:
  1. Reordered original data matrix
  2. Indices giving the new order of points
"""
function _reorder_within_clusters(data::Matrix{Float64}, 
                                data_quantiles::Matrix{Float64}, 
                                clusters::Vector{Int},
                                params_quantile::Vector{Float64})
    n_clusters = maximum(clusters)
    n_cols = size(data, 2)
    
    # Initialize output arrays
    reordered_data = copy(data)
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
            
            # Find the point closest to params_quantile
            distances_to_params = [sum((cluster_points[:, i] .- params_quantile).^2) for i in 1:n_points]
            current_idx = argmin(distances_to_params)
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
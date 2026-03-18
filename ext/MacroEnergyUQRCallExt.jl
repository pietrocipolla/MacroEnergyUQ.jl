module MacroEnergyUQRCallExt

using MacroEnergyUQ, RCall

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
function MacroEnergyUQ.transport_ot(source_marginal::Vector{Float64}, target_marginal::Vector{Float64}, cost_matrix::Matrix{Float64})
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
function MacroEnergyUQ.transport_ranking(marginal::Vector{Float64}, cost_matrix::Matrix{Float64})
    return MacroEnergyUQ.transport_ot(marginal, marginal, cost_matrix)
end

"""
    ot_indices(x::Matrix, y::Matrix, M::Int; kwargs...)

Wrapper around the R package gsaot to compute Optimal Transport-based sensitivity indices.

# Arguments
- `x`: Input matrix where each row is an observation and each column is a variable. Values can be numeric, strings, or categorical.
       For continuous inputs (Float64), data is partitioned into M sets. For discrete inputs, partitioning is data-driven.
- `y`: Output matrix where each row is an observation and each column is an output variable. Must be numeric.
- `M`: Number of partitions for continuous inputs.
- `cost`: Cost function for the Optimal Transport problem. Default is "L2" (squared Euclidean metric).
         Can be either "L2" or a custom function that takes y as input and returns a cost matrix.
- `discrete_out`: If true, creates a histogram of realizations to use as weights. Useful for discrete or mixed outputs
                 with large number of realizations. Default is false.
- `solver`: Solver for the Optimal Transport problem. Options:
           - "sinkhorn" (default): Sinkhorn's solver
           - "sinkhorn_log": Sinkhorn's solver in log scale
           - "transport": Non-regularized OT problem solver
- `solver_optns`: Optional dictionary of solver options.
- `scaling`: Whether to scale the cost matrix. Default is true.
- `boot`: Whether to perform bootstrapping of OT indices. Default is false.
- `stratified_boot`: If true (default), creates partitions first then performs stratified bootstrapping.
                    If false, resamples dataset before creating partitions.
- `R`: Number of bootstrap replicas. Default is nothing (no bootstrapping).
- `parallel`: Type of parallel operation for bootstrapping. Options: "no" (default), "multicore", "snow".
- `ncpus`: Number of CPUs for parallel operation. Default is 1.
- `conf`: Confidence level for bootstrap confidence intervals (0-1). Default is 0.95.
- `type`: Method for computing confidence intervals. Default is "norm".
         Other options: "basic", "stud", "perc", "bca".

# Returns
A dictionary containing:
- Target sensitivity indices
- Confidence intervals (if bootstrapping is enabled)
"""
function MacroEnergyUQ.ot_indices(x::Matrix, y::Matrix, M::Int; 
                   cost::String="L2",
                   discrete_out::Bool=false,
                   solver::String="sinkhorn",
                   solver_optns=nothing,
                   scaling::Bool=true,
                   boot::Bool=false,
                   stratified_boot::Bool=true,
                   R=nothing,
                   parallel::String="no",
                   ncpus::Int=1,
                   conf::Float64=0.95,
                   type::String="norm")

    # Check gsaot package is available in R
    R"""
    if (!require("gsaot", quietly = TRUE)) {
        stop("R package 'gsaot' is not installed. Please install it with: install.packages('gsaot')")
    }
    """
    
    # Convert Julia arrays to R objects
    @rput x
    @rput y
    @rput M
    @rput cost
    @rput discrete_out
    @rput solver
    @rput solver_optns
    @rput scaling
    @rput boot
    @rput stratified_boot
    @rput R
    @rput parallel
    @rput ncpus
    @rput conf
    @rput type

    R"""
    # Ensure X is a matrix
    x <- as.data.frame(t(x))
    
    # Compute target sensitivity indices using optimal transport
    result <- gsaot::ot_indices(x = x, y = y, M = M, 
                               cost = cost,
                               discrete_out = discrete_out,
                               solver = solver,
                               solver_optns = solver_optns,
                               scaling = scaling,
                               boot = boot,
                               stratified_boot = stratified_boot,
                               R = R,
                               parallel = parallel,
                               ncpus = ncpus,
                               conf = conf,
                               type = type)
    """

    # Get results back to Julia
    @rget result

    return result
end

"""
    ot_indices_wb(x::Matrix, y::Matrix, M::Int; kwargs...)

Wrapper around the R package gsaot to compute Wasserstein-Bures sensitivity indices.

# Arguments
- `x`: Input matrix where each row is an observation and each column is a variable. Values can be numeric, strings, or categorical.
       For continuous inputs (Float64), data is partitioned into M sets. For discrete inputs, partitioning is data-driven.
- `y`: Output matrix where each row is an observation and each column is an output variable. Must be numeric.
- `M`: Number of partitions for continuous inputs.
- `boot`: Whether to perform bootstrapping of OT indices. Default is false.
- `R`: Number of bootstrap replicas. Default is nothing (no bootstrapping).
- `parallel`: Type of parallel operation for bootstrapping. Options: "no" (default), "multicore", "snow".
- `ncpus`: Number of CPUs for parallel operation. Default is 1.
- `conf`: Confidence level for bootstrap confidence intervals (0-1). Default is 0.95.
- `type`: Method for computing confidence intervals. Default is "norm".
         Other options: "basic", "stud", "perc", "bca".

# Returns
A dictionary containing:
- Target sensitivity indices
- Confidence intervals (if bootstrapping is enabled)
"""
function MacroEnergyUQ.ot_indices_wb(x::Matrix, y::Matrix, M::Int; 
                       boot::Bool=false,
                       R=nothing,
                       parallel::String="no",
                       ncpus::Int=1,
                       conf::Float64=0.95,
                       type::String="norm")
    
    # Convert Julia arrays to R objects
    @rput x
    @rput y
    @rput M
    @rput boot
    @rput R
    @rput parallel
    @rput ncpus
    @rput conf
    @rput type

    R"""
    # Ensure X is a matrix
    x <- as.data.frame(t(x))
    
    # Compute target sensitivity indices using optimal transport
    result <- gsaot::ot_indices_wb(x = x, y = y, M = M, 
                                   boot = boot,
                                   R = R,
                                   parallel = parallel,
                                   ncpus = ncpus,
                                   conf = conf,
                                   type = type)
    """

    # Get results back to Julia
    @rget result

    return result
end

"""
    ot_indices_1d(x::Matrix, y::Array, M::Int; kwargs...)

Wrapper around the R package gsaot to compute Optimal Transport sensitivity indices for univariate outputs.

# Arguments
- `x`: Input matrix where each row is an observation and each column is a variable. Values can be numeric, strings, or categorical.
       For continuous inputs (Float64), data is partitioned into M sets. For discrete inputs, partitioning is data-driven.
- `y`: Output array where each row is an observation and each column is an output variable. Must be numeric.
- `M`: Number of partitions for continuous inputs.
- `boot`: Whether to perform bootstrapping of OT indices. Default is false.
- `R`: Number of bootstrap replicas. Default is nothing (no bootstrapping).
- `parallel`: Type of parallel operation for bootstrapping. Options: "no" (default), "multicore", "snow".
- `ncpus`: Number of CPUs for parallel operation. Default is 1.
- `conf`: Confidence level for bootstrap confidence intervals (0-1). Default is 0.95.
- `type`: Method for computing confidence intervals. Default is "norm".
         Other options: "basic", "stud", "perc", "bca".

# Returns
A dictionary containing:
- Target sensitivity indices
- Confidence intervals (if bootstrapping is enabled)
"""
function MacroEnergyUQ.ot_indices_1d(x::Matrix, y::Array, M::Int;
                       p::Float64 = 2, 
                       boot::Bool=false,
                       R=nothing,
                       parallel::String="no",
                       ncpus::Int=1,
                       conf::Float64=0.95,
                       type::String="norm")
    
    # Convert Julia arrays to R objects
    @rput x
    @rput y
    @rput M
    @rput p
    @rput boot
    @rput R
    @rput parallel
    @rput ncpus
    @rput conf
    @rput type

    R"""
    # Ensure X is a matrix
    x <- as.data.frame(t(x))
    
    # Compute target sensitivity indices using optimal transport
    result <- gsaot::ot_indices_1d(x = x, y = y, M = M, p = p,
                                   boot = boot,
                                   R = R,
                                   parallel = parallel,
                                   ncpus = ncpus,
                                   conf = conf,
                                   type = type)
    """

    # Get results back to Julia
    @rget result

    return result
end


"""
    irrelevance_threshold(y, M; dummy_optns=nothing, cost="L2", discrete_out=false, 
                         solver="sinkhorn", solver_optns=nothing, scaling=true)

Wrapper around the R package gsaot to calculate irrelevance threshold using dummy variable 
for Optimal Transport sensitivity indices.

# Arguments
- `y`: An array or a matrix containing the output values.
- `M`: A scalar representing the number of partitions for continuous inputs.
- `dummy_optns`: (default nothing) A list containing the options on the distribution of the dummy variable.
- `cost`: (default "L2") A string or function defining the cost function of the Optimal Transport problem. 
         It should be "L2" or a function taking as input y and returning a cost matrix. 
         If cost="L2", uses the squared Euclidean metric.
- `discrete_out`: (default false) Logical, by default the output sample in y are equally weighted. 
                 If discrete_out=true, the function tries to create a histogram of the realizations 
                 and to use the histogram as weights. Useful for discrete or mixed outputs with large 
                 number of realizations. Reduces the dimension of the cost matrix.
- `solver`: Solver for the Optimal Transport problem. Options:
           - "1d": one-dimensional analytic solution
           - "wasserstein-bures": Wasserstein-Bures solution
           - "sinkhorn" (default): Sinkhorn's solver
           - "sinkhorn_log": Sinkhorn's solver in log scale
           - "transport": non-regularized OT problem solver using transport::transport()
- `solver_optns`: (optional) A list containing the options for the Optimal Transport solver.
- `scaling`: (default true) Whether or not to scale the cost matrix.

# Returns
A dictionary containing the irrelevance threshold results.
"""
function MacroEnergyUQ.irrelevance_threshold(y::Matrix, M::Int; 
                              dummy_optns=nothing,
                              cost::String="L2",
                              discrete_out::Bool=false,
                              solver::String="sinkhorn",
                              solver_optns=nothing,
                              scaling::Bool=true)
    
    # Convert Julia arrays to R objects
    @rput y
    @rput M
    @rput dummy_optns
    @rput cost
    @rput discrete_out
    @rput solver
    @rput solver_optns
    @rput scaling
    
    R"""
    # Calculate irrelevance threshold using dummy variable
    result <- gsaot::irrelevance_threshold(
        y = y,
        M = M,
        dummy_optns = dummy_optns,
        cost = cost,
        discrete_out = discrete_out,
        solver = solver,
        solver_optns = solver_optns,
        scaling = scaling
    )
    """
    
    # Get results back to Julia
    @rget result
    
    return result
end

end

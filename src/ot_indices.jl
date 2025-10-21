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
function ot_indices(x::Matrix, y::Matrix, M::Int; 
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
function ot_indices_wb(x::Matrix, y::Matrix, M::Int; 
                    cost::String="L2",
                    discrete_out::Bool=false,
                    solver::String="sinkhorn",
                    solver_optns=nothing,
                    scaling::Bool=true,
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
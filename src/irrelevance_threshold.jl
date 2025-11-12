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
function irrelevance_threshold(y::Matrix, M::Int; 
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

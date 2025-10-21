"""
    run_mc(model_factory::Function, data::Matrix{Float64}, 
           params::AbstractVector{String},
           optimizer;
           clusters::Vector{Int} = ones(Int, size(data, 2)),
           kwargs...)

Run Monte Carlo simulations using a function that creates a JuMP model. The function creates a new model for each thread, updates the model parameters for each sample in `data`, solves the model, and collects results. The simulations are parallelized across threads.

# Arguments
- `model_factory`: A function that returns a new JuMP model. The function should accept an optimizer and additional keyword arguments
- `data`: Matrix where each column represents a parameter set for one simulation
- `params`: Vector of variable names corresponding to the uncertain parameters in the model (must match the number of rows in `data`)
- `optimizer`: The optimizer to use (can be a type like `Gurobi.Optimizer` or a function)
- `clusters`: (Optional) A vector indicating the cluster assignment for each sample in `data`
- `kwargs...`: Additional keyword arguments that will be passed to `model_factory`


# Returns
A NamedTuple containing:
- `status`: Vector of optimization status for each simulation
- `solve_time`: Vector of solve times for each simulation
- `outputs`: Matrix where each row corresponds to a solution vector

# Example
```julia
using JuMP, MacroEnergyUQ

# Create a simple model
model = Model()
@variable(model, x >= 0)
@variable(model, y >= 0)
@objective(model, Min, 0.1*x + 0.6*y)
@constraint(model, x + y >= 1)

# Generate sample points
data = QuasiMonteCarlo.sample(100, 2, SobolSample())

# Run Monte Carlo simulations
results = run_mc(model, data, params=["x", "y"])
status = results.status      # Vector of solver status
times = results.solve_time   # Vector of computation times
Y = results.outputs         # Matrix of solutions
```
"""
function run_mc(model_factory::Function, data::Matrix{Float64}, 
                params::AbstractVector{String},
                optimizer;
                clusters::Vector{Int} = ones(Int, size(data, 2)),
                kwargs...)
    # Validate inputs
    if length(params) != size(data, 1)
        error("Number of parameters ($(length(params))) does not match the number of rows in data ($(size(data, 1))).")
    end

    if length(clusters) != size(data, 2)
        error("Length of clusters ($(length(clusters))) must match the number of samples in data ($(size(data, 2))).")
    end
    
    n_samples = size(data, 2)
    n_threads = maximum(clusters)
    
    # Prepare storage for results for each thread
    thread_results = [Vector{Dict{Symbol, Any}}() for _ in 1:n_threads]
    
    # Create base models for each thread
    thread_models = [model_factory(optimizer; kwargs...) for _ in 1:n_threads]
    
    # Run simulations in parallel
    Threads.@threads for cluster in 1:n_threads
        # Get indices for this cluster
        cluster_indices = findall(==(cluster), clusters)
        
        # Get the model and results for this thread
        thread_model = thread_models[cluster]
        cluster_results = Dict{Symbol, Any}[]
        
        # Process each sample in this cluster
        for idx in cluster_indices
            # Update parameters
            for (param_idx, param_name) in enumerate(params)
                set_objective_coefficient(thread_model, 
                                          variable_by_name(thread_model, param_name), 
                                          data[param_idx, idx])
            end
            
            # Solve the model and store results
            time_start = time()
            optimize!(thread_model)
            solve_time = time() - time_start
            
            # Store results in thread-local array
            push!(cluster_results, Dict{Symbol, Any}(
                :index => idx,
                :status => termination_status(thread_model),
                :objective => objective_value(thread_model),
                :solution => value.(all_variables(thread_model)),
                :solve_time => solve_time
            ))

            thread_results[cluster] = cluster_results
        end
    end
    
    # Initialize arrays for results
    status = Vector{Any}(undef, n_samples)
    solve_time = Vector{Float64}(undef, n_samples)
    
    # Get the size of the solution vector from the first result
    first_solution = thread_results[1][1][:solution]
    outputs = Matrix{Float64}(undef, n_samples, length(first_solution))
    
    # Combine results from all threads
    for thread_result in thread_results
        for cluster_result in thread_result
            idx = cluster_result[:index]
            status[idx] = cluster_result[:status]
            solve_time[idx] = cluster_result[:solve_time]
            outputs[idx, :] = cluster_result[:solution]
        end
    end
    
    return (status = status, solve_time = solve_time, outputs = outputs)
end
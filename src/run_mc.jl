"""
    run_mc(model::Model, data::Matrix{Float64}; 
           params::AbstractVector{Symbol} = Symbol[], 
           n_threads::Int = Threads.nthreads(),
           optimizer_factory = nothing)

Run Monte Carlo simulations on a JuMP model using the provided dataset. The function updates the model parameters in-place for each sample in `data`, solves the model, and collects results. The simulations are parallelized across `n_threads` threads.

# Arguments
- `model`: A JuMP model
- `data`: Matrix where each column represents a parameter set for one simulation
- `params`: Vector of variable names (as symbols) corresponding to the uncertain parameters in the model (must match the number of rows in `data`)

# Returns
- `results`: A vector of dictionaries, each containing:
  - `:status`: Optimization status
  - `:objective`: Objective value (if solved)
  - `:solution`: Solution values (if solved)
  - `:solve_time`: Time taken to solve the model

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
```
"""
function run_mc(model::Model, data::Matrix{Float64}, 
                params::AbstractVector{String};
                clusters::Vector{Int} = ones(Int, size(data, 2)))
    
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
    thread_models = [copy(model) for _ in 1:n_threads]
    
    # Run simulations in parallel
    Threads.@threads for cluster in 1:n_threads
        # Get indices for this cluster
        cluster_indices = findall(==(cluster), clusters)
        
        # Get the model and results for this thread
        thread_model = thread_models[cluster]
        cluster_results = thread_results[cluster]
        
        # Process each sample in this cluster
        for idx in cluster_indices
            # Update parameters
            for (param_idx, param_name) in enumerate(params)
                set_objective_coefficient(thread_model, 
                                          variable_by_name(thread_model, param_name), 
                                          ordered_data[param_idx, idx])
            end
            
            # Solve the model and store results
            time_start = time()
            optimize!(thread_model)
            solve_time = time() - time_start
            
            # Check solution status
            assert_is_solved_and_feasible(model)
            
            # Store results in thread-local array
            push!(cluster_results, Dict{Symbol, Any}(
                :index => idx,
                :status => termination_status(thread_model),
                :objective => objective_value(thread_model),
                :solution => value.(all_variables(thread_model)),
                :solve_time => solve_time
            ))
        end
    end
    
    # Combine results from all threads
    results = Vector{Dict{Symbol, Any}}(undef, n_samples)
    for cluster_results in thread_results
        for result in cluster_results
            results[result[:index]] = result
        end
        # Remove the :index field from the final results
        for result in values(results)
            delete!(result, :index)
        end
    end
    
    return results
end
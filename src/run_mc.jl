function run_cluster(cluster::Int, 
                model_factory::Function, 
                data::Matrix{Float64}, 
                params::AbstractVector{String},
                optimizer;
                params_type::Symbol = :objective,
                clusters::Vector{Int} = ones(Int, size(data, 2)),
                extract::Function = m -> value.(all_variables(m)),
                kwargs...)

  # Get indices for this cluster
  cluster_indices = findall(==(cluster), clusters)
        
  # Get the model and results for this thread
  model = model_factory(optimizer; kwargs...)
  cluster_results = Dict{Symbol, Any}[]
        
  # Process each sample in this cluster
  for idx in cluster_indices
    # Update parameters
    for (param_idx, param_name) in enumerate(params)
      if params_type == :objective
        set_objective_coefficient(model, 
                                  variable_by_name(model, param_name), 
                                  data[param_idx, idx])
      else
        set_parameter_value(variable_by_name(model, param_name), 
                            data[param_idx, idx])
      end
    end
            
    # Solve and extract (error handling delegated to extract function)
    time_start = time()
    optimize!(model)
    solve_time = time() - time_start
    
    push!(cluster_results, Dict{Symbol, Any}(
          :cluster => cluster,
          :index => idx,
          :status => termination_status(model),
          :solution => extract(model),
          :solve_time => solve_time
          ))
  end

  return cluster_results
end

"""
    run_mc(model_factory::Function, data::Matrix{Float64}, 
           params::AbstractVector{String}, optimizer;
           params_type::Symbol = :objective,
           clusters::Vector{Int} = ones(Int, size(data, 2)),
           extract::Function = m -> value.(all_variables(m)),
           kwargs...)

Run parallel Monte Carlo simulations using a JuMP model factory.

`run_mc` creates a new JuMP model for each thread (or cluster), updates
the model parameters for each Monte Carlo sample, solves the model, and
collects user-specified outputs. Simulations are parallelized across threads.

# Arguments
- `model_factory`: A function that returns a new JuMP model.
  The function should accept an optimizer and additional keyword arguments.

- `data`: A matrix where each **column** represents a sample of uncertain
  parameters, and each **row** corresponds to a parameter in `params`.

- `params`: Vector of variable names (as `String`s) corresponding to the uncertain
  parameters in the model. Must match the number of rows in `data`.

- `optimizer`: The optimizer to use (e.g., `Gurobi.Optimizer`, `Ipopt.Optimizer`).

- `params_type`: A symbol representing how to find the parameters in the model. If
  `params_type == :objective`, the strings in `params` are assumed to be the name of the
  variable corresponding to the uncertain coefficient, and `set_objective_coefficient` is used.
  If `params_type == :parameter`, the strings in `params` are assumed to be the name of the
  uncertain parameters themselves, and `set_parameter_value` is used.

- `clusters`: (Optional) A vector of integers indicating the cluster assignment
  for each sample in `data`. Samples assigned to the same cluster are solved
  sequentially on the same thread. Defaults to all ones (single cluster).

- `extract`: A function that extracts outputs from the solved model.
  It should accept a JuMP model and return a vector (or any indexable object)
  of numerical results.
  
  By default, `extract = m -> value.(all_variables(m))` (returns all variable values).

  Example custom extractor:
  ```julia
  extract = m -> [
      value(m[:x]),                   # single variable
      value(m[:y]),                   # another variable
      value(0.1*m[:x] + 0.6*m[:y])   # custom expression
  ]

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
                params_type::Symbol = :objective,
                clusters::Vector{Int} = ones(Int, size(data, 2)),
                extract::Function = m -> value.(all_variables(m)),
                distributed::Bool = false,
                kwargs...)
    # Validate inputs
    if length(params) != size(data, 1)
        error("Number of parameters ($(length(params))) does not match the number of rows in data ($(size(data, 1))).")
    end

    if length(clusters) != size(data, 2)
        error("Length of clusters ($(length(clusters))) must match the number of samples in data ($(size(data, 2))).")
    end

    if params_type != :objective && params_type != :parameter
      throw(ArgumentError("Unsupported reference to parameters: $params_type. Currently, only :objective and :parameter are supported."))
    end
    
    n_samples = size(data, 2)
    n_processes = maximum(clusters)

    mc_results_raw = [Vector{Dict{Symbol, Any}}() for _ in 1:n_processes]
    
    # Run distributed simulations
    mc_results_raw = pmap(cluster -> run_cluster(cluster, model_factory, data, params, optimizer;
                                                params_type=params_type,
                                                clusters=clusters,
                                                extract=extract,
                                                kwargs...), 
                          1:n_processes; 
                          distributed=distributed)
    
    # Initialize arrays for results
    status = Vector{Any}(undef, n_samples)
    solve_time = Vector{Float64}(undef, n_samples)
    
    # Get the size of the solution vector from the first result
    first_solution = mc_results_raw[1][1][:solution]
    outputs = Matrix{Float64}(undef, n_samples, length(first_solution))
    
    # Combine results from all threads, using the stored index to maintain input order
    for mc_result_raw in mc_results_raw
        for cluster_result in mc_result_raw
            idx = cluster_result[:index]
            status[idx] = cluster_result[:status]
            solve_time[idx] = cluster_result[:solve_time]
            outputs[idx, :] = cluster_result[:solution]
        end
    end

    return (status = status, solve_time = solve_time, outputs = outputs)
end
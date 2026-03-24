# Helper function to update parameters in a model
function update_parameters!(model::Model, 
                           params::AbstractVector{String}, 
                           data::Matrix{Float64}, 
                           param_idx_col::Int, 
                           params_type::Symbol)
  for (param_idx, param_name) in enumerate(params)
    if params_type == :objective
      set_objective_coefficient(model,
                                variable_by_name(model, param_name),
                                data[param_idx, param_idx_col])
    else
      set_parameter_value(variable_by_name(model, param_name),
                          data[param_idx, param_idx_col])
    end
  end
end

# Helper function to save or store results
function save_result!(cluster_results::Vector{Dict{Symbol, Any}}, 
                     cluster::Int, 
                     idx::Int, status, 
                     solve_time::Float64,
                     solution::Vector, 
                     write_outputs::Bool, 
                     output_dir::String)
  result_dict = Dict{Symbol, Any}(
    :cluster => cluster,
    :index => idx,
    :status => status,
    :solve_time => solve_time
  )
  
  if write_outputs
    output_file = joinpath(output_dir, "sample_$(idx).csv")
    CSV.write(output_file, Tables.table(solution'), writeheader=false)
    result_dict[:solution_file] = output_file
  else
    result_dict[:solution] = solution
  end
  
  push!(cluster_results, result_dict)
end

# Process cluster for standard JuMP models (without context)
function process_cluster_samples(model::Model, 
                                cluster_indices::Vector{Int},
                                cluster::Int, 
                                data::Matrix{Float64},
                                params::AbstractVector{String}, 
                                params_type::Symbol,
                                extract::Function, 
                                write_outputs::Bool, 
                                output_dir::String)
  cluster_results = Dict{Symbol, Any}[]
  
  for idx in cluster_indices
    update_parameters!(model, params, data, idx, params_type)
    
    time_start = time()
    optimize!(model)
    solve_time = time() - time_start
    
    solution = extract(model)
    save_result!(cluster_results, cluster, idx, termination_status(model), 
                solve_time, solution, write_outputs, output_dir)
  end
  
  return cluster_results
end

# Process cluster for standard JuMP models (with context)
function process_cluster_samples(model::Model, 
                                context::NamedTuple,
                                cluster_indices::Vector{Int},
                                cluster::Int, 
                                data::Matrix{Float64},
                                params::AbstractVector{String}, 
                                params_type::Symbol,
                                extract::Function, 
                                write_outputs::Bool, 
                                output_dir::String)
  cluster_results = Dict{Symbol, Any}[]
  
  for idx in cluster_indices
    update_parameters!(model, params, data, idx, params_type)
    
    time_start = time()
    optimize!(model)
    solve_time = time() - time_start
    
    # Add index to context for this sample
    ctx_with_index = merge(context, (index=idx,))
    solution = extract(model; ctx=ctx_with_index)
    save_result!(cluster_results, cluster, idx, termination_status(model), 
                solve_time, solution, write_outputs, output_dir)
  end
  
  return cluster_results
end

# Process cluster for Benders decomposition (without context)
function process_cluster_samples(components::NamedTuple, 
                                cluster_indices::Vector{Int},
                                cluster::Int, 
                                data::Matrix{Float64},
                                params::AbstractVector{String}, 
                                params_type::Symbol,
                                extract::Function, 
                                write_outputs::Bool, 
                                output_dir::String)
  cluster_results = Dict{Symbol, Any}[]
  planning_problem = components.planning_problem
  benders_settings = get(components, :settings, Dict())
  
  # Check if this has context - if so, dispatch to context version
  if haskey(components, :context)
    components_context = components.context
    return process_cluster_samples(components, components_context, 
                                  cluster_indices, cluster, data, params, 
                                  params_type, extract, write_outputs, output_dir)
  end
  
  for idx in cluster_indices
    update_parameters!(planning_problem, params, data, idx, params_type)
    
    time_start = time()
    results = MacroEnergySolvers.benders(
      components.planning_problem,
      components.subproblems,
      components.linking_variables,
      benders_settings
    )
    solve_time = time() - time_start
    
    solution = extract(results)
    save_result!(cluster_results, cluster, idx, "BENDERS", 
                solve_time, solution, write_outputs, output_dir)
  end
  
  return cluster_results
end

# Process cluster for Benders decomposition (with context)
function process_cluster_samples(components::NamedTuple, 
                                context::NamedTuple,
                                cluster_indices::Vector{Int},
                                cluster::Int, 
                                data::Matrix{Float64},
                                params::AbstractVector{String}, 
                                params_type::Symbol,
                                extract::Function, 
                                write_outputs::Bool, 
                                output_dir::String)
  cluster_results = Dict{Symbol, Any}[]
  planning_problem = components.planning_problem
  benders_settings = get(components, :settings, Dict())
  
  for idx in cluster_indices
    update_parameters!(planning_problem, params, data, idx, params_type)
    
    time_start = time()
    results = MacroEnergySolvers.benders(
      components.planning_problem,
      components.subproblems,
      components.linking_variables,
      benders_settings
    )
    solve_time = time() - time_start
    
    # Add index to context for this sample
    ctx_with_index = merge(context, (index=idx,))
    solution = extract(results; ctx=ctx_with_index)
    save_result!(cluster_results, cluster, idx, "BENDERS", 
                solve_time, solution, write_outputs, output_dir)
  end
  
  return cluster_results
end

# Main cluster processing function with multiple dispatch
function run_cluster(cluster::Int, 
                    model_factory::Function, 
                    data::Matrix{Float64}, 
                    params::AbstractVector{String},
                    optimizer;
                    params_type::Symbol = :objective,
                    clusters::Vector{Int} = ones(Int, size(data, 2)),
                    extract::Function = m -> value.(all_variables(m)),
                    write_outputs::Bool = false,
                    cluster_start_delay_s::Float64 = 0.0,
                    output_dir::String = ".",
                    kwargs...)

  cluster_indices = findall(==(cluster), clusters)

  if cluster_start_delay_s > 0.0
    sleep((cluster - 1) * cluster_start_delay_s)
  end

  components = model_factory(optimizer; kwargs...)
  
  # Determine if this is a standard JuMP model (has :model key) or Benders components
  if haskey(components, :model)
    # Standard JuMP model: extract model and context
    model = components.model
    context = get(components, :context, nothing)
    
    if context !== nothing
      return process_cluster_samples(model, context, cluster_indices, cluster, 
                                    data, params, params_type, extract, 
                                    write_outputs, output_dir)
    else
      return process_cluster_samples(model, cluster_indices, cluster, 
                                    data, params, params_type, extract, 
                                    write_outputs, output_dir)
    end
  else
    # Benders components: dispatch with or without context
    context = get(components, :context, nothing)
    
    if context !== nothing
      return process_cluster_samples(components, context, cluster_indices, cluster, 
                                    data, params, params_type, extract, 
                                    write_outputs, output_dir)
    else
      return process_cluster_samples(components, cluster_indices, cluster, 
                                    data, params, params_type, extract, 
                                    write_outputs, output_dir)
    end
  end
end

"""
    run_mc(model_factory::Function, data::Matrix{Float64}, 
           params::AbstractVector{String}, optimizer;
           params_type::Symbol = :objective,
           clusters::Vector{Int} = ones(Int, size(data, 2)),
           extract::Function = m -> value.(all_variables(m)),
           distributed::Bool = false,
           write_outputs::Bool = false,
           output_dir::String = "mc_outputs",
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

- `extract`: A function that extracts outputs from the solved model/results.
  It can accept either one argument (the solved model/results) or two arguments 
  (with a `ctx` keyword argument containing context data).
  
  The context is automatically provided if `model_factory` returns a NamedTuple
  with a `context` field. The context can contain fields like `index` (sample index)
  and any other data the extract function needs (e.g., the full system for post-processing).
  
  By default, `extract = m -> value.(all_variables(m))` (returns all variable values).

  Example custom extractors:
  ```julia
  # Simple extractor without context
  extract = m -> [
      value(m[:x]),                   # single variable
      value(m[:y]),                   # another variable
      value(0.1*m[:x] + 0.6*m[:y])   # custom expression
  ]
  
  # Extractor with context (receives ctx keyword)
  extract = (results; ctx=nothing) -> [
      value(results.planning_sol.objective_value),
      # ...use ctx.index for output directory
      # ...use ctx.system for post-processing
  ]
  ```

- `distributed`: (Optional) Whether to use distributed computing with `pmap`.
  Defaults to `false` (uses multi-threading).

- `write_outputs`: (Optional) If `true`, solution vectors are written to individual
  CSV files in `output_dir` instead of being stored in memory. This is useful for
  large models where storing all solutions in memory is impractical. 
  Defaults to `false`.

- `output_dir`: (Optional) Directory path where solution files will be written when
  `write_outputs=true`. Each solution is saved as `sample_<idx>.csv`. The directory
  is created automatically if it doesn't exist. Defaults to `"mc_outputs"`.

- `kwargs...`: Additional keyword arguments that will be passed to `model_factory`


# Returns
A NamedTuple containing:
- `status`: Vector of optimization status for each simulation
- `solve_time`: Vector of solve times for each simulation
- `outputs`: Matrix where each row corresponds to a solution vector (only when `write_outputs=false`)
- `output_dir`: Path to directory containing solution files (only when `write_outputs=true`)

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
                write_outputs::Bool = false,
                cluster_start_delay_s::Float64 = 0.0,
                output_dir::String = "mc_outputs",
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
    
    # Create output directory if writing outputs
    if write_outputs
        mkpath(output_dir)
    end
    
    n_samples = size(data, 2)
    n_processes = maximum(clusters)

    mc_results_raw = [Vector{Dict{Symbol, Any}}() for _ in 1:n_processes]
    
    # Run distributed simulations
    mc_results_raw = pmap(cluster -> run_cluster(cluster, model_factory, data, params, optimizer;
                                                params_type=params_type,
                                                clusters=clusters,
                                                extract=extract,
                                                write_outputs=write_outputs,
                                                cluster_start_delay_s=cluster_start_delay_s,
                                                output_dir=output_dir,
                                                kwargs...), 
                          1:n_processes; 
                          distributed=distributed)
    
    # Initialize arrays for results
    status = Vector{Any}(undef, n_samples)
    solve_time = Vector{Float64}(undef, n_samples)
    
    if write_outputs
        # Return only metadata when outputs are written to disk
        for mc_result_raw in mc_results_raw
            for cluster_result in mc_result_raw
                idx = cluster_result[:index]
                status[idx] = cluster_result[:status]
                solve_time[idx] = cluster_result[:solve_time]
            end
        end
        
        return (status = status, solve_time = solve_time, output_dir = output_dir)
    else
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
end
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

# Helper functions for Benders cuts
"""
    name_cuts(planning_problem::Model, next_cut_id::Int) -> Int

Assign deterministic names (`Cut_<id>`) to unnamed constraints.
Returns the next available cut id.
"""
function name_cuts(planning_problem::Model, next_cut_id::Int)
  for con in all_constraints(planning_problem, include_variable_in_set_constraints=false)
    if isempty(name(con))
      set_name(con, "Cut_" * string(next_cut_id))
      next_cut_id += 1
    end
  end
  return next_cut_id
end

# Slack computed according to the tutorial https://jump.dev/JuMP.jl/dev/tutorials/linear/lp_sensitivity/#Constraint-sensitivity
function cut_slack(con::JuMP.ConstraintRef)
  slack = normalized_rhs(con) - value(con)
  # Numerical safeguard: tiny negative slacks can appear due to tolerances.
  if slack < 0.0
    return abs(slack)
  end
  return slack
end

function cut_constraints(planning_problem::Model)
  filter(con -> startswith(name(con), "Cut_"),
         all_constraints(planning_problem, include_variable_in_set_constraints=false))
end

function is_feasibility_cut(planning_problem::Model, con::JuMP.ConstraintRef)
  if !haskey(object_dictionary(planning_problem), :vTHETA)
    return false
  end

  theta_ref = planning_problem[:vTHETA]
  theta_vars = theta_ref isa AbstractArray ? vec(theta_ref) : [theta_ref]

  for vtheta in theta_vars
    coeff = normalized_coefficient(con, vtheta)
    if abs(coeff) > 1e-12
      return false
    end
  end

  return true
end

"""
    manage_cuts(planning_problem, constraints_to_keep, benders_settings)

Apply cut storage policy and remove constraints that should not be retained.

Supported `CutSetting` values:
- `"store_all"`: keep all constraints.
- `"store_lc"`: keep only constraints present after the first Benders iteration.
- `"store_first_n"`: keep only cuts `Cut_1` to `Cut_n` (`n = benders_settings[:num_cuts]`, default 10000).
- `"store_none"`: keep only the base constraints (those present before cut naming starts).
- `"store_binding"`: keep base constraints and only cuts that are binding (or recent, via aging).

Optional settings for `"store_binding"`:
- `BindingTol` (default `1e-7`): slack tolerance below which a cut is considered binding.
- `MinCutAge` (default `0`): minimum number of iterations to keep a cut regardless of binding status.
- `MaxInactiveAge` (default `0`): number of consecutive non-binding iterations a cut can survive.
"""
function manage_cuts(planning_problem::Model,
                     constraints_to_keep::Vector{String},
                     benders_settings,
                     cut_age::Union{Nothing, Dict{String, Int}},
                     cut_inactive_age::Union{Nothing, Dict{String, Int}})

  cut_setting = get(benders_settings, :CutSetting, "store_all")

  if cut_setting == "store_all"
    return name.(cut_constraints(planning_problem))
  elseif cut_setting == "store_lc"
    # First invocation stores the current constraint set. Subsequent invocations
    # remove newly generated cuts and keep only the stored set.
    cuts_already_saved = !isempty(constraints_to_keep)
    if cuts_already_saved
      forget_cuts!(planning_problem, constraints_to_keep)
    else
      return name.(cut_constraints(planning_problem))
    end
  elseif cut_setting == "store_first_n"
    num_cuts = get(benders_settings, :NumCuts, 10000)
    @debug "Managing Benders cuts with 'store_first_n'" num_cuts
    for con in cut_constraints(planning_problem)
      con_name = name(con)
      cut_idx_match = match(r"^Cut_(\d+)$", con_name)
      if cut_idx_match !== nothing
        cut_idx = parse(Int, cut_idx_match.captures[1])
        if cut_idx <= num_cuts && !(con_name in constraints_to_keep)
          push!(constraints_to_keep, con_name)
        end
      end
    end
    forget_cuts!(planning_problem, constraints_to_keep)
    return constraints_to_keep
  elseif cut_setting == "store_none"
    forget_cuts!(planning_problem, constraints_to_keep)
  elseif cut_setting == "store_binding"
    # Read parameters for binding-based storage
    binding_tol = get(benders_settings, :BindingTol, 1e-7)
    min_cut_age = get(benders_settings, :MinCutAge, 1)
    max_inactive_age = get(benders_settings, :MaxInactiveAge, 2)

    keep = Set(constraints_to_keep)
    active_cut_names = Set{String}()

    for con in cut_constraints(planning_problem)
      con_name = name(con)
      push!(active_cut_names, con_name)
      current_age = get(cut_age, con_name, 0) + 1
      cut_age[con_name] = current_age

      # Feasibility cuts are structural safeguards and should never be pruned
      # by binding/inactive logic, otherwise infeasible subproblems can recur.
      if is_feasibility_cut(planning_problem, con)
        cut_inactive_age[con_name] = 0
        push!(keep, con_name)
        continue
      end

      if current_age <= min_cut_age
        cut_inactive_age[con_name] = 0
        push!(keep, con_name)
        continue
      end

      slack = cut_slack(con)
      if slack <= binding_tol
        cut_inactive_age[con_name] = 0
        push!(keep, con_name)
      else
        inactive_age = get(cut_inactive_age, con_name, 0) + 1
        cut_inactive_age[con_name] = inactive_age
        if inactive_age <= max_inactive_age
          push!(keep, con_name)
        else
          delete!(keep, con_name)
        end
      end
    end

    # Remove stale entries from the aging dictionaries.
    for cut_name in collect(keys(cut_age))
      if !(cut_name in active_cut_names)
        delete!(cut_age, cut_name)
        delete!(cut_inactive_age, cut_name)
      end
    end

    # Drop references to cuts that no longer exist in the planning model.
    for cut_name in collect(keep)
      if startswith(cut_name, "Cut_") && !(cut_name in active_cut_names)
        delete!(keep, cut_name)
      end
    end

    constraints_to_keep = collect(keep)
    forget_cuts!(planning_problem, constraints_to_keep)
  else
    error("Invalid cut storage setting: $cut_setting. Must be one of 'store_all', 'store_lc', 'store_first_n', 'store_none', or 'store_binding'.")
  end

  return constraints_to_keep
end

function forget_cuts!(planning_problem::Model, constraints_to_keep::Vector{String})
  keep = Set(constraints_to_keep)
  to_delete = JuMP.ConstraintRef[]

  for con in cut_constraints(planning_problem)
    if !(name(con) in keep)
      push!(to_delete, con)
    end
  end

  @info "Removing $(length(to_delete)) cuts from the planning problem"

  for con in to_delete
    delete(planning_problem, con)
  end
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

  counter = 1
  planning_constraints = name.(cut_constraints(planning_problem))
  cut_age = Dict{String, Int}()
  cut_inactive_age = Dict{String, Int}()
  
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
    save_result!(cluster_results, cluster, idx, results.termination_status, 
                solve_time, solution, write_outputs, output_dir)

    ### Cut Management ###

    cut_setting = get(benders_settings, :CutSetting, "store_all")

    if cut_setting == "store_binding"
      # set_name(...) is a model modification and can invalidate primal values.
      # For binding-based filtering, evaluate slack before assigning new cut names.
      planning_constraints = manage_cuts(planning_problem, planning_constraints, benders_settings,
                                         cut_age, cut_inactive_age)
      counter = name_cuts(planning_problem, counter)
    else
      # Label cuts by iteration
      counter = name_cuts(planning_problem, counter)

      # Manage cut storage
      planning_constraints = manage_cuts(planning_problem, planning_constraints, benders_settings,
                                         cut_age, cut_inactive_age)
    end

    ### End Cut Management ###
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


  counter = 1
  planning_constraints = name.(cut_constraints(planning_problem))
  cut_age = Dict{String, Int}()
  cut_inactive_age = Dict{String, Int}()
  
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
    save_result!(cluster_results, cluster, idx, results.termination_status, 
                solve_time, solution, write_outputs, output_dir)

    ### Cut Management ###

    cut_setting = get(benders_settings, :CutSetting, "store_all")

    if cut_setting == "store_binding"
      # set_name(...) is a model modification and can invalidate primal values.
      # For binding-based filtering, evaluate slack before assigning new cut names.
      planning_constraints = manage_cuts(planning_problem, planning_constraints, benders_settings,
                                         cut_age, cut_inactive_age)
      counter = name_cuts(planning_problem, counter)
    else
      # Label cuts by iteration
      counter = name_cuts(planning_problem, counter)

      # Manage cut storage
      planning_constraints = manage_cuts(planning_problem, planning_constraints, benders_settings,
                                         cut_age, cut_inactive_age)
    end

    ### End Cut Management ###
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
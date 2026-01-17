"""
    run_benders(master_factory, subproblem_factory, distributed, write_outputs, master_optimizer; kwargs...)
    Run Benders decomposition to find optimal solution. 


"""

function run_benders(master_factory::Function, # Function to build master problem ----- must have original objective function labeled as :eObj and theta variables labeled as :vTHETA
    subproblem_factory::Function, # Function to build subproblems
    distributed::Bool = false,
    write_outputs::Bool = false,
    master_optimizer;
    subproblem_optimizer;
    workers::Int = 1,
    inputs_decomp::Dict, # Decomposed inputs for subproblems - time domain split into rep periods equal to number of subproblems
    outputs_dir::String = "benders_outputs",
    extractor_subproblem::Function = m -> value.(all_variables(m)), # Function to extract subproblem solutions
	extractor_master::Function = m -> value.(all_variables(m)); # Function to extract master solution
    kwargs...
)
    if distributed
        using Distributed
        n_workers = nworkers()
        if n_workers < workers
            addprocs(workers - n_workers)
        end
    end

    # Build models
    planning_problem, planning_variables  = master_factory(master_optimizer; kwargs...)
    subproblems, planning_variables_sub = subproblem_factory(subproblem_optimizer, inputs_decomp; kwargs...)

    # Initialize Benders inputs
    benders_inputs = Dict{String, Any}()

    # Update benders inputs
    benders_inputs["planning_problem"] = planning_problem
    benders_inputs["planning_variables"] = planning_variables
    benders_inputs["subproblems"] = subproblems
    benders_inputs["planning_variables_sub"] = planning_variables_sub

    # Run Benders decomposition
    benders_result = benders(
        benders_inputs,
        extractor_subproblem, # Function to extract subproblem solutions
        extractor_master, # Function to extract master solution
        inputs,
        setup,
        kwargs...
    )

    return benders_result

end

#=============================================================================
    Run Benders with Multiple MGA Iterations


=============================================================================#
"""
    run_benders_mga(master_factory, subproblem_factory, inputs, setup; kwargs...)

Run Benders decomposition to find optimal solution, then perform multiple MGA iterations.

# Arguments
- `master_factory::Function`: Function `(optimizer; kwargs...) -> (model, variables)` to build master problem
- `subproblem_factory::Function`: Function `(optimizer; kwargs...) -> (subproblems, variables_sub)` to build subproblems
- `inputs::Dict`: Problem-specific input data
- `setup::Dict`: Algorithm configuration parameters

# Keyword Arguments
## Optimizers
- `master_optimizer`: Optimizer for the master problem
- `subproblem_optimizer`: Optimizer for subproblems

## MGA Configuration
- `n_iterations::Int`: Number of MGA iterations to run (default: 1)
- `objective_factory::Function`: Function `(model, iteration) -> expr` to build iteration-specific objective

## Model Building (optional)
- `extractor_subproblem::Function`: Function to extract subproblem solutions
- `extractor_master::Function`: Function to extract master solutions

## MGA Customization (optional)
- `setup_mga_problem!::Function`: Function `(model, setup) -> nothing` to set up MGA constraints
- `name_cuts!::Function`: Function `(model, counter) -> counter` to name existing cuts
- `extract_master_solution::Function`: Function `(model, vars, inputs, metadata...) -> NamedTuple` 
- `extract_results::Function`: Function `(master_sol, subop_sol) -> Any` to format results
- `update_cuts!::Function`: Function to add Benders cuts
- `solve_subproblems::Function`: Function to solve subproblems

## Distributed Computing (optional)
- `distributed::Bool`: Whether to use distributed computing (default: false)
- `workers::Int`: Number of workers (default: 1)

## Output (optional)
- `write_outputs::Bool`: Whether to write outputs (default: false)
- `outputs_dir::String`: Directory for outputs (default: "benders_outputs")

# Returns
- `opt_result`: Results from the initial optimization
- `mga_results`: Array of results for each MGA iteration
- `timing_df`: DataFrame with timing statistics
"""
function run_benders_mga(
    master_factory::Function,
    subproblem_factory::Function,
    inputs::Dict,
    setup::Dict;
    # Required arguments
    master_optimizer,
    subproblem_optimizer,
    n_iterations::Int = 1,
    objective_factory::Function,
    # Optional model extractors
    extractor_subproblem::Function = m -> value.(all_variables(m)),
    extractor_master::Function = m -> value.(all_variables(m)),
    # Optional MGA/investment sensitivity customization
    setup_mga_problem!::Function = (model, setup) -> nothing,
    name_cuts!::Function = default_name_cuts!,
    extract_results::Function = (master_sol, subop_sol) -> (master_sol, subop_sol),
    # Distributed computing
    distributed::Bool = false,
    workers::Int = 1,
    # Output options
    write_outputs::Bool = false,
    outputs_dir::String = "benders_outputs",
    inputs_decomp::Dict = Dict(),
    kwargs...
)
    println("=" ^ 60)
    println("Phase 1: Running initial Benders optimization")
    println("=" ^ 60)
    
    # Check for distributed computing
    if distributed
        using Distributed
        n_workers = nworkers()
        if n_workers < workers
            addprocs(workers - n_workers)
        end
    end

    # Build models using factory functions
    planning_problem, planning_variables = master_factory(master_optimizer; kwargs...)
    subproblems, planning_variables_sub = subproblem_factory(subproblem_optimizer, inputs_decomp; kwargs...)

    # Initialize Benders inputs
    benders_inputs = Dict{String, Any}()
    benders_inputs["planning_problem"] = planning_problem
    benders_inputs["planning_variables"] = planning_variables
    benders_inputs["subproblems"] = subproblems
    benders_inputs["planning_variables_sub"] = planning_variables_sub

    # Run initial Benders decomposition to find optimal solution
    opt_result = benders(
        benders_inputs,
        extractor_master,
        extractor_subproblem,
        inputs,
        setup;
        kwargs...
    )
    
    # Extract models, subproblems, and variables from optimization result
    planning_problem = opt_result.benders_inputs["planning_problem"]
    planning_variables = opt_result.benders_inputs["planning_variables"]
    subproblems = opt_result.benders_inputs["subproblems"]
    planning_variables_sub = opt_result.benders_inputs["planning_variables_sub"]


    nsubs = length(subproblems)
    
    println("\n" * "=" ^ 60)
    println("Phase 2: Running MGA iterations ($n_iterations iterations)")
    println("=" ^ 60)
    
    # Calculate MGA budget from optimal solution
    setup["MGABudget"] = opt_result.UB_hist[end] * (1 + setup["ModelingtoGenerateAlternativeSlack"])
    println("MGA Budget: ", setup["MGABudget"])

    # Set up MGA master problem (user-defined or no-op)
    setup_mga_problem!(planning_problem, setup)

    # Name and track existing cuts
    cut_counter = 0
    cut_counter = name_cuts!(planning_problem, cut_counter)
    opt_cuts = name.(all_constraints(planning_problem, include_variable_in_set_constraints = false))

    sumtime_df = DataFrame(
        :MGA_it => 0, 
        :Iterations => length(opt_result.UB_hist), 
        :Iteration_Time => opt_result.cpu_time[end]
    )

    retain_master_cuts = get(setup, "ModelingToGenerateAlternativeRetainBendersCuts", 2)
    println("Cut Retention Setting: ", retain_master_cuts)
    setup["BD_Stab_Method"] = "off"
    
    mga_results = Vector{Any}(undef, n_iterations)

    for iteration in 1:n_iterations
        println("\n--- MGA Iteration $iteration of $n_iterations ---")
        
        # Handle cut retention strategy
        _apply_cut_retention!(planning_problem, opt_cuts, retain_master_cuts, setup, nsubs)
        
        # Set iteration-specific objective using factory function
        obj_expr = objective_factory(planning_problem, iteration)
        @objective(planning_problem, Min, obj_expr)

        # Run MGA cutting plane algorithm
        @time mga_result = mga_cutting_plane(
                benders_inputs,
                setup, 
                inputs, 
                iteration;
                extractor_master,
                extractor_subproblem
            )
        
        # Update planning problem reference (in case it was modified)
        planning_problem = mga_result.EP_master
        master_sol_final = mga_result.master_sol
        subop_sol_mga = mga_result.subop_sol
        
        # Extract and store results using user-defined function
        mga_results[iteration] = extract_results(master_sol_final, subop_sol_mga)
        
        # Log results
        _log_iteration_results(master_sol_final, subop_sol_mga)
    
        time_df = DataFrame(
            :MGA_it => iteration, 
            :Iterations => length(mga_result.TrueSystemCost_hist), 
            :Iteration_Time => mga_result.cpu_time[end]
        )
        append!(sumtime_df, time_df)
    end
    
    println("\n" * "=" ^ 60)
    println("Completed: 1 optimal + $n_iterations MGA iterations")
    println("=" ^ 60)
    
    return opt_result, mga_results, sumtime_df
end


function run_benders_investment_sensitivity(
    master_factory::Function,
    subproblem_factory::Function,
    inputs::Dict,
    setup::Dict;
    # Required arguments
    master_optimizer,
    subproblem_optimizer,
    investment_sensitivity_values::Vector{Dict{String, Any}},
    objective_factory::Function,
    # Optional model extractors
    extractor_subproblem::Function = m -> value.(all_variables(m)),
    extractor_master::Function = m -> value.(all_variables(m)),
    # Distributed computing
    distributed::Bool = false,
    workers::Int = 1,
    # Output options
    write_outputs::Bool = false,
    outputs_dir::String = "benders_outputs",
    inputs_decomp::Dict = Dict(),
    kwargs...
)
    println("=" ^ 60)
    println("Running Investment Sensitivity Analysis")
    println("=" ^ 60)
    
    # Check for distributed computing
    if distributed
        using Distributed
        n_workers = nworkers()
        if n_workers < workers
            addprocs(workers - n_workers)
        end
    end

    # Build models using factory functions
    planning_problem, planning_variables = master_factory(master_optimizer; kwargs...)
    subproblems, planning_variables_sub = subproblem_factory(subproblem_optimizer, inputs_decomp; kwargs...)

    # Initialize Benders inputs
    benders_inputs = Dict{String, Any}()
    benders_inputs["planning_problem"] = planning_problem
    benders_inputs["planning_variables"] = planning_variables
    benders_inputs["subproblems"] = subproblems
    benders_inputs["planning_variables_sub"] = planning_variables_sub

    # Store results for each sensitivity value
    sensitivity_results = Vector{Any}(undef, length(investment_sensitivity_values))
    
    for (idx, value) in enumerate(investment_sensitivity_values)
        println("\n--- Investment Sensitivity Iteration $idx ---")
        
        update_objective = objective_factory(planning_problem, value)
        @objective(planning_problem, Min, update_objective)

        # Run Benders decomposition for current sensitivity value
        sensitivity_result = benders(
            benders_inputs,
            extractor_subproblem,
            extractor_master,
            inputs,
            setup;
            kwargs...
        )
        sensitivity_results[idx] = sensitivity_result
    end


    println("\n" * "=" ^ 60)
    println("Completed Investment Sensitivity Analysis")
    println("=" ^ 60)

    return sensitivity_results
end

#=============================================================================
    Helper Functions for Benders Running Functions
=============================================================================#

"""
Default function to name cuts in the master problem.
"""
function default_name_cuts!(model::Model, counter::Int)
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        if name(con) == ""
            set_name(con, "BendersCut" * string(counter))
            counter += 1
        end
    end
    return counter
end

"""
Apply cut retention strategy based on setup configuration.
"""
function _apply_cut_retention!(model::Model, opt_cuts::Vector{String}, retain_method::Int, setup::Dict, nsubs::Int)
    if retain_method == 1
        # Keep all cuts
    elseif retain_method == 2
        forget_cuts_master!(model, opt_cuts)
    elseif retain_method == 3
        recent_cuts = retain_recent_cuts(model, opt_cuts, setup["MaxCuts"])
        forget_cuts_master!(model, recent_cuts)
    elseif retain_method == 5
        sp_cuts = retain_fixed_spcuts_early(model, opt_cuts, setup["MaxCuts"], nsubs)
        forget_cuts_master!(model, sp_cuts)
    else
        println("No cut-retention method specified, defaulting to optimal cuts only")
        forget_cuts_master!(model, opt_cuts)
    end
end

"""
Default function to extract master problem solution.
Override with custom extractor for model-specific fields.
"""
function default_extract_master_solution(model::Model, master_vars::Vector{String}, inputs::Dict, id::Int, iteration::Int, mga_it::Int)
    return (
        inv_cost = value(model[:eObj]),
        values = Dict([s => value.(variable_by_name(model, s)) for s in master_vars]),
        id = id,
        iteration = iteration,
        mga_it = mga_it
    )
end

"""
Log iteration results (can be overridden or extended).
"""
function _log_iteration_results(master_sol::NamedTuple, subop_sol::Dict)
    println("Investment cost: ", master_sol.inv_cost)
    if haskey(master_sol, :zone_inv_cost)
        println("Zone investment cost: ", master_sol.zone_inv_cost)
    end
    println("Operational costs: ", [subop_sol[w].op_cost for w in keys(subop_sol)])
    if haskey(first(values(subop_sol)), :zone_cost)
        println("Zone operational costs: ", [subop_sol[w].zone_cost for w in keys(subop_sol)])
    end
end

function retain_recent_cuts(model::Model, master_cons::Vector{String}, num_cuts::Int)
    cut_names = Vector{String}()
    opt_names = Vector{String}()
    struc_names = Vector{String}()
    
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        con_name = name(con)
        if con_name == "" || occursin("BendersCut", con_name)
            split_name = split(con_name, "_")
            if length(split_name) >= 2
                mga_it = tryparse(Int, split_name[2])
                if mga_it !== nothing && mga_it == 0
                    push!(opt_names, con_name)
                else
                    push!(cut_names, con_name)
                end
            end
        else
            push!(struc_names, con_name)
        end
    end
    
    opt_count = length(opt_names)
    total_cuts = length(cut_names)
    start_idx = max(1, total_cuts - num_cuts - opt_count)
    
    retained = [opt_names; cut_names[start_idx:end]]
    return [struc_names; retained]
end


function retain_fixed_spcuts_early(model::Model, master_cons::Vector{String}, num_cuts::Int, nsubs::Int)
    sp_cuts = [Vector{String}() for _ in 1:nsubs]
    struc_names = Vector{String}()
    
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        con_name = name(con)
        if occursin("BendersCut", con_name)
            split_name = split(con_name, "_")
            if length(split_name) >= 4
                num_str = split(split_name[4], "[")
                sp_idx = tryparse(Int, num_str[1])
                if sp_idx !== nothing && 1 <= sp_idx <= nsubs
                    push!(sp_cuts[sp_idx], con_name)
                end
            end
        else
            push!(struc_names, con_name)
        end
    end
    
    cut_names = Vector{String}()
    for i in 1:nsubs
        n_cuts = min(length(sp_cuts[i]), num_cuts)
        append!(cut_names, sp_cuts[i][1:n_cuts])
    end
    
    return [struc_names; cut_names]
end

function forget_cuts_master!(model::Model, retained_cons::Vector{String})
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        if !(name(con) in retained_cons)
            delete(model, con)
        end
    end
end

"""
    solve_mga_master_problem(model, master_vars, inputs, id, iteration, mga_it; kwargs...)

Solve the MGA master problem and extract solution.

# Arguments
- `model::Model`: The master problem JuMP model
- `master_vars::Vector{String}`: Names of master variables to extract
- `inputs::Dict`: Problem inputs
- `id::Int`: Solution identifier
- `iteration::Int`: Current iteration number
- `mga_it::Int`: MGA iteration number

# Keyword Arguments
- `extract_master_solution::Function`: Custom solution extractor function
- `validate_solution::Function`: Optional function to validate solution (returns bool)
- `crossover_fallback::Bool`: Whether to retry with crossover on negative values

# Returns
- `NamedTuple`: Master problem solution
"""
function solve_mga_master_problem(
    model::Model,
    master_vars::Vector{String}, 
    inputs::Dict, 
    id::Int, 
    iteration::Int, 
    mga_it::Int;
    extract_master_solution::Function = default_extract_master_solution,
    validate_solution::Function = (m) -> true,
    crossover_fallback::Bool = true
)
    iteration += 1
    println("Solving MGA master problem...")
    optimize!(model)
    
    # Check if solution needs crossover fallback
    needs_crossover = !validate_solution(model)
    
    if needs_crossover && crossover_fallback
        println("***Resolving master problem with Crossover=1***")
        set_attribute(model, "Crossover", 1)
        optimize!(model)
        if has_values(model)
            master_sol = extract_master_solution(model, master_vars, inputs, id, iteration, mga_it)
            set_attribute(model, "Crossover", 0)
        else
            error("Master problem failed to solve even with crossover enabled")
        end
    else
        master_sol = extract_master_solution(model, master_vars, inputs, id, iteration, mga_it)
    end
    
    return master_sol
end



"""
Check if the MGA budget constraint is satisfied.
"""
function _check_budget_convergence(cost::Float64, setup::Dict)
    budget = setup["MGABudget"]
    relax = get(setup, "RelaxBudget", 0.0)
    
    if relax > 0 && isapprox(cost, budget, rtol=relax)
        return true
    end
    return cost <= budget
end

"""
Log timing statistics for the MGA iteration.
"""
function _log_timing_stats(master_times::Vector{Float64}, sub_times::Vector{Float64})
    if !isempty(master_times) && !isempty(sub_times)
        master_avg = mean(master_times)
        subop_avg = mean(sub_times)
        ms_ratio = master_avg / subop_avg
        println("MGA iteration finished")
        println("Average Master Time = $master_avg")
        println("Average Subop Time = $subop_avg")
        println("Master/Subop Ratio = $ms_ratio")
    end
end

"""
    update_master_problem_multi_cuts_mga!(model, subop_sol, master_sol, master_vars_sub, mga_it, k)

Default function to add Benders cuts to the master problem for MGA.
Override this for model-specific cut generation.
"""
function update_master_problem_multi_cuts_mga!(model::Model, subop_sol::Dict, master_sol::NamedTuple, master_vars_sub::Dict, mga_it::Int, k::Int)
    W = keys(subop_sol)
    
    for w in W
        cut_name = "BendersCut_$(mga_it)_$(k)_$w"
        @constraint(model,
            subop_sol[w].theta_coeff * model[:vTHETA][w] >= 
            subop_sol[w].op_cost + sum(
                subop_sol[w].lambda[i] * (variable_by_name(model, master_vars_sub[w][i]) - master_sol.values[master_vars_sub[w][i]]) 
                for i in 1:length(master_vars_sub[w])
            ),
            base_name = cut_name
        )
    end
end


#=============================================================================
    Helper Functions for Model-Specific Adapters
=============================================================================#

"""
    make_linear_objective_factory(coefficients::Dict{Symbol, Array})

Create an objective factory function from a dictionary of variable symbols and their coefficient arrays.

# Arguments
- `coefficients`: Dict mapping variable symbols to coefficient arrays indexed by iteration

# Returns
- Function `(model, iteration) -> AffExpr` suitable for `objective_factory` argument

# Example
```julia
# For a model with :vCap[tech, zone] and :vLine[line] variables
coefficients = Dict(
    :vCap => rand(n_techs, n_zones, n_iterations),
    :vLine => rand(n_lines, n_iterations)
)
obj_factory = make_linear_objective_factory(coefficients)
```
"""
function make_linear_objective_factory(coefficients::Dict{Symbol, <:AbstractArray})
    return function(model::Model, iteration::Int)
        expr = AffExpr(0.0)
        for (var_sym, coef_array) in coefficients
            var = model[var_sym]
            dims = ndims(coef_array)
            if dims == 2  # (index, iteration)
                n = size(coef_array, 1)
                add_to_expression!(expr, sum(coef_array[i, iteration] * var[i] for i in 1:n))
            elseif dims == 3  # (index1, index2, iteration)
                n1, n2 = size(coef_array, 1), size(coef_array, 2)
                add_to_expression!(expr, sum(coef_array[i, j, iteration] * var[i, j] for i in 1:n1, j in 1:n2))
            else
                error("Unsupported coefficient array dimension: $dims")
            end
        end
        return expr
    end
end

"""
    make_capacity_validation_function(capacity_vars::Vector{Symbol})

Create a validation function that checks for negative capacity values.

# Arguments
- `capacity_vars`: Vector of symbols for capacity variables to check

# Returns
- Function `(model) -> Bool` suitable for `validate_solution` argument
"""
function make_capacity_validation_function(capacity_vars::Vector{Symbol})
    return function(model::Model)
        for var_sym in capacity_vars
            if haskey(model, var_sym)
                if any(value.(model[var_sym]) .< 0)
                    return false
                end
            end
        end
        return true
    end
end

"""
    make_solution_extractor(inv_cost_expr::Symbol, additional_fields::Dict{Symbol, Function})

Create a master solution extractor with customizable fields.

# Arguments
- `inv_cost_expr`: Symbol for the investment cost expression in the model
- `additional_fields`: Dict of field names to extractor functions `(model, inputs) -> value`

# Returns  
- Function suitable for `extract_master_solution` argument
"""
function make_solution_extractor(inv_cost_expr::Symbol, additional_fields::Dict{Symbol, Function} = Dict{Symbol, Function}())
    return function(model::Model, master_vars::Vector{String}, inputs::Dict, id::Int, iteration::Int, mga_it::Int)
        base_sol = (
            inv_cost = value(model[inv_cost_expr]),
            values = Dict([s => value.(variable_by_name(model, s)) for s in master_vars]),
            id = id,
            iteration = iteration,
            mga_it = mga_it
        )
        
        # Add additional fields
        extra = Dict{Symbol, Any}()
        for (field_name, extractor) in additional_fields
            extra[field_name] = extractor(model, inputs)
        end
        
        if isempty(extra)
            return base_sol
        else
            return merge(base_sol, NamedTuple(extra))
        end
    end
end

"""
    setup_budget_constraint!(model::Model, setup::Dict, theta_var::Symbol = :vTHETA, obj_expr::Symbol = :eObj)

Add a budget constraint for MGA to the master problem.

# Arguments
- `model`: The JuMP model
- `setup`: Setup dictionary containing "MGABudget" key
- `theta_var`: Symbol for theta variables (default: :vTHETA)
- `obj_expr`: Symbol for objective expression (default: :eObj)
"""
function setup_budget_constraint!(model::Model, setup::Dict; theta_var::Symbol = :vTHETA, obj_expr::Symbol = :eObj)
    @constraint(model, cMGABudget, model[obj_expr] + sum(model[theta_var]) == setup["MGABudget"])
end



function name_cuts!(EP_master::Model, counter::Int64)
    for con in all_constraints(EP_master,include_variable_in_set_constraints=false)
        if name(con) == ""
            set_name(con,"BendersCut"*string(counter))
        end
        counter+=1
    end 
    return counter
end

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
function default_extract_master_solution(model::Model, master_vars::Vector{String})
    return (
        inv_cost = value(model[:eObj]),
        values = Dict([s => value.(variable_by_name(model, s)) for s in master_vars])
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


function setup_mga_master_problem!(EP_master::Model,setup::Dict)
 

    @constraint(EP_master,cMGABudget, EP_master[:eObj] + sum(EP_master[:vTHETA]) == setup["MGABudget"])

end

#=============================================================================
    Helper Functions for Benders Running Functions
=============================================================================#


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

function default_objective_factory!(planning_problem::Model, mga_variables::Vector{String}, vector::Vector{Float64}; method::Int = 0, kwargs...)
    if length(vector) != length(mga_variables)
        error("Length of objective vector must match number of MGA variables")
    end
    @objective(planning_problem, Min, sum(vector[i] * variable_by_name(planning_problem, mga_variables[i]) for i in 1:length(mga_variables)))
    
    return
end

function default_vector_generator(mga_variables::Vector{String}, iterations::Int; method::Int = 0, kwargs...)
    vars = length(mga_variables)
    coeffs = zeros(iterations, vars)
    if method == 0
        coeffs = rand(Float64, (iterations, vars))
    elseif method == 1
        coeffs = rand(-1:1, (iterations, vars))
    elseif method == 2
        coeffs_1 = rand(Float64, (ceil(Int64,iterations*0.25), vars))
        coeffs_2 = rand(-1:1, (floor(Int64,iterations*0.75), vars))
        coeffs = vcat(coeffs_1, coeffs_2)
    end
    return coeffs
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

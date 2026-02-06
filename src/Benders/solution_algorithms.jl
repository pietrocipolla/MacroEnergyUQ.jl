function benders(
	benders_inputs::Dict{String, Any},
    inputs::Dict,
    setup::Dict;
    extractor_subproblem::Function = m -> value.(all_variables(m)), # Function to extract subproblem solutions
	extractor_master::Function = m -> value.(all_variables(m)),
    check_negative_capacities::Function = m -> default_check_negative_capacities(m), # Function to extract master solution;
    kwargs...)


    # Set algorithm parameters
    MaxIter = setup["BD_MaxIter"]
    ConvTol = setup["BD_ConvTol"]
	MaxCpuTime = setup["BD_MaxCpuTime"]
	γ = setup["BD_StabParam"];
	stab_method = setup["BD_Stab_Method"];
    integer_investment = setup["IntegerInvestments"]
    integer_routine_flag = false

	#### Retrieve models
	planning_problem = benders_inputs["planning_problem"]
	planning_variables = benders_inputs["planning_variables"]
	subproblems = benders_inputs["subproblems"]
	planning_variables_sub = benders_inputs["planning_variables_sub"]

    if integer_investment == 1 && stab_method != "off"
		all_planning_variables = all_variables(planning_problem);
		integer_variables = all_planning_variables[is_integer.(all_planning_variables)];
		binary_variables = all_planning_variables[is_binary.(all_planning_variables)];
		unset_integer.(integer_variables)
		unset_binary.(binary_variables)
		integer_routine_flag = true;
	end

    solver_start_time = time()

    #### Initialize UB and LB
	planning_sol = MacroEnergyUQ.solve_planning_problem(planning_problem,planning_variables, extractor_master, check_negative_capacities);
	subop_sol = Dict()

    UB = Inf;
    LB = planning_sol.LB;

    LB_hist = Float64[];
    UB_hist = Float64[];
    cpu_time = Float64[];
	feasibility_hist = Float64[];

	planning_sol_best = deepcopy(planning_sol);

    #### Run Benders iterations
    for k = 0:MaxIter
		
		start_subop_sol = time();

        subop_sol = solve_dist_subproblems(subproblems,planning_sol,inputs, extractor_subproblem);
        
		cpu_subop_sol = time()-start_subop_sol;
		println("Solving the subproblems required $cpu_subop_sol seconds")

		UBnew = sum((subop_sol[w].theta_coeff==0 ? Inf : subop_sol[w].op_cost) for w in keys(subop_sol))+planning_sol.inv_cost;
		if UBnew < UB
			planning_sol_best = deepcopy(planning_sol);
			UB = UBnew;
		end

		print("Updating the planning problem....")
		time_start_update = time()

		update_planning_problem_multi_cuts!(planning_problem,subop_sol,planning_sol,planning_variables_sub)
		
		time_planning_update = time()-time_start_update
		println("done (it took $time_planning_update s).")

		start_planning_sol = time()
		unst_planning_sol = solve_planning_problem(planning_problem,planning_variables, extractor_master, check_negative_capacities);
		cpu_planning_sol = time()-start_planning_sol;
		println("Solving the planning problem required $cpu_planning_sol seconds")

		LB = max(LB,unst_planning_sol.LB);
		
		append!(LB_hist,LB)
        append!(UB_hist,UB)
		append!(feasibility_hist,sum(subop_sol[w].feasibility_slack for w in keys(subop_sol)))
        append!(cpu_time,time()-solver_start_time)

		if any(subop_sol[w].theta_coeff==0 for w in keys(subop_sol))
			println("***k = ", k,"      LB = ", LB,"     UB = ", UB,"       Gap = ", (UB-LB)/abs(LB),"       CPU Time = ",cpu_time[end])
		else
			println("k = ", k,"      LB = ", LB,"     UB = ", UB,"       Gap = ", (UB-LB)/abs(LB),"       CPU Time = ",cpu_time[end])
		end

        if (UB-LB)/abs(LB) <= ConvTol
			if integer_routine_flag
				println("*** Switching on integer constraints *** ")
				UB = Inf;
				set_integer.(integer_variables)
				set_binary.(binary_variables)
				planning_sol = solve_planning_problem(planning_problem,planning_variables, extractor_master, check_negative_capacities);
				LB = planning_sol.LB;
				planning_sol_best = deepcopy(planning_sol);
				integer_routine_flag = false;
			else
				break
			end
		elseif (cpu_time[end] >= MaxCpuTime)|| (k == MaxIter)
			break
		elseif UB==Inf
			planning_sol = deepcopy(unst_planning_sol);
		else
			if stab_method == "int_level_set"
				start_stab_method = time()
				if  integer_investment==1 && integer_routine_flag==false
					unset_integer.(integer_variables)
					unset_binary.(binary_variables)
					for v in integer_variables
						fix(v,unst_planning_sol.values[name(v)];force=true)
					end
					for v in binary_variables
						fix(v,unst_planning_sol.values[name(v)];force=true)
					end
                    println("Solving the interior level set problem with γ = $γ")
					planning_sol = solve_int_level_set_problem(planning_problem,planning_variables,unst_planning_sol,LB,UB,γ,inputs, extractor_master);
					unfix.(integer_variables)
					unfix.(binary_variables)
					set_integer.(integer_variables)
					set_binary.(binary_variables)
					set_lower_bound.(integer_variables,0.0)
					set_lower_bound.(binary_variables,0.0)
				else
                    println("Solving the interior level set problem with γ = $γ")
					planning_sol = solve_int_level_set_problem(planning_problem,planning_variables,unst_planning_sol,LB,UB,γ,inputs, extractor_master);
				end
				cpu_stab_method = time()-start_stab_method;
				println("Solving the interior level set problem required $cpu_stab_method seconds")
			else
				planning_sol = deepcopy(unst_planning_sol);
			end

		end

    end

	### Update Benders inputs
	benders_inputs["planning_problem"] = planning_problem
	benders_inputs["planning_variables"] = planning_variables 
	benders_inputs["subproblems"] = subproblems
	benders_inputs["planning_variables_sub"] = planning_variables_sub

	return (
        benders_inputs = benders_inputs,
        planning_sol = planning_sol_best,
        operational_sol = subop_sol,
        LB_hist = LB_hist,
        UB_hist = UB_hist,
        cpu_time = cpu_time,
        feasibility_hist = feasibility_hist
    )
end

"""
    mga_cutting_plane(model, master_vars, subproblems, master_vars_sub, setup, inputs, mga_it; kwargs...)

Run the MGA cutting plane algorithm for a single MGA iteration.

# Arguments
- `benders_inputs::Dict{String, Any}`: Benders algorithm inputs - contains master problem, subproblems, and variable lists
- `setup::Dict`: Algorithm configuration
- `inputs::Dict`: Problem inputs
- `mga_it::Int`: MGA iteration number

# Keyword Arguments
- `extract_master_solution::Function`: Custom solution extractor
- `update_cuts!::Function`: Function to add cuts to master problem
- `solve_subproblems::Function`: Function to solve subproblems

# Returns
- Named tuple with master problem, solutions, cost histories, and timing
"""
function mga_cutting_plane(
    benders_inputs::Dict{String, Any},
    setup::Dict, 
    inputs::Dict, 
    mga_it::Int;
    extract_master_solution::Function = default_extract_master_solution,
	extractor_subproblem::Function = m -> value.(all_variables(m))
)
    # Initialize timing
    cpu_time = [0.0]
    solver_start_time = time()
    id = 1
    indicator = 0

    # Algorithm parameters
    MaxIter = setup["BD_MaxIter"]
    MaxCpuTime = setup["BD_MaxCpuTime"]

    TrueSystemCost = Inf
    ApproxSystemCost = setup["MGABudget"]

    ApproxSystemCost_hist = [ApproxSystemCost]
    TrueSystemCost_hist = [TrueSystemCost]
    
    master_sol_final = (inv_cost = 0.0, values = Dict{String,Any}())
    subop_sol = Dict()
    
    master_times = Float64[]
    sub_times = Float64[]

	#### Retrieve models
	planning_problem = benders_inputs["planning_problem"]
	planning_variables = benders_inputs["planning_variables"]
	subproblems = benders_inputs["subproblems"]
	planning_variables_sub = benders_inputs["planning_variables_sub"]

    for k = 1:MaxIter
        # Solve master problem
        start_master_sol = time()
        master_sol = solve_mga_master_problem(
            planning_problem, planning_variables, inputs, id, k, mga_it;
            extract_master_solution = extract_master_solution
        )
        cpu_master_sol = time() - start_master_sol
        println("Solving the master problem required $cpu_master_sol seconds")

        # Solve subproblems
        start_subop_sol = time()
        subop_sol = solve_subproblems(subproblems, master_sol, inputs, extractor_subproblem)
        cpu_subop_sol = time() - start_subop_sol
        push!(sub_times, cpu_subop_sol)
        println("Solving the subproblems required $cpu_subop_sol seconds")

        # Update best solution
        TrueSystemCost_new = sum(subop_sol[w].op_cost for w in keys(subop_sol)) + master_sol.inv_cost
        if TrueSystemCost_new <= TrueSystemCost
            TrueSystemCost = TrueSystemCost_new
            master_sol_final = deepcopy(master_sol)
        end

        push!(ApproxSystemCost_hist, ApproxSystemCost)
        push!(TrueSystemCost_hist, TrueSystemCost)
        push!(cpu_time, time() - solver_start_time)

        println("k = $k, ApproxSystemCost = $ApproxSystemCost, TrueSystemCost = $TrueSystemCost, ",
                "TrueSystemCost_new = $TrueSystemCost_new, ",
                "MGABudget Violation = $((TrueSystemCost_new - setup["MGABudget"]) / abs(setup["MGABudget"])), ",
                "CPU Time = $(cpu_time[end])")

        # Check convergence
        budget_satisfied = _check_budget_convergence(TrueSystemCost_new, setup)
        
        if budget_satisfied
            if indicator == 0
                println("Rerunning with crossover on")
                set_attribute(model, "Crossover", 1)
                TrueSystemCost = Inf
                indicator = 1
            else
                set_attribute(model, "Crossover", 0)
                _log_timing_stats(master_times, sub_times)
                return (
                    EP_master = model,
                    master_sol = master_sol_final,
                    subop_sol = subop_sol,
                    ApproxSystemCost_hist = ApproxSystemCost_hist,
                    TrueSystemCost_hist = TrueSystemCost_hist,
                    cpu_time = cpu_time
                )
            end
        elseif cpu_time[end] >= MaxCpuTime
            return (
                EP_master = model,
                master_sol = master_sol_final,
                subop_sol = subop_sol,
                ApproxSystemCost_hist = ApproxSystemCost_hist,
                TrueSystemCost_hist = TrueSystemCost_hist,
                cpu_time = cpu_time
            )
        else
            # Update master problem with cuts
            print("Updating the master problem....")
            time_start_update = time()
            
            update_cuts!(model, subop_sol, master_sol, master_vars_sub, mga_it, k)
            
            time_master_update = time() - time_start_update
            println("done (it took $time_master_update s).")
            push!(master_times, cpu_master_sol + time_master_update)
        end
    end

	### Update Benders inputs
	benders_inputs["planning_problem"] = planning_problem
	benders_inputs["planning_variables"] = planning_variables 
	benders_inputs["subproblems"] = subproblems
	benders_inputs["planning_variables_sub"] = planning_variables_sub

    # Return after max iterations
    return (
        benders_inputs = benders_inputs,
        master_sol = master_sol_final,
        subop_sol = subop_sol,
        ApproxSystemCost_hist = ApproxSystemCost_hist,
        TrueSystemCost_hist = TrueSystemCost_hist,
        cpu_time = cpu_time
    )
end

function update_planning_problem_multi_cuts!(EP::Model,subop_sol::Dict,planning_sol::NamedTuple,planning_variables_sub::Dict)
    
	W = keys(subop_sol);
	
    @constraint(EP,[w in W],subop_sol[w].theta_coeff*EP[:vTHETA][w] >= subop_sol[w].op_cost + sum(subop_sol[w].lambda[i]*(variable_by_name(EP,planning_variables_sub[w][i]) - planning_sol.values[planning_variables_sub[w][i]]) for i in 1:length(planning_variables_sub[w])));

       
end
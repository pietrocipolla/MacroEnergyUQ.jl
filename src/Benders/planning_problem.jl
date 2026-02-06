


function solve_planning_problem(EP::Model,planning_variables::Vector{String},extractor_planning_problem::Function, check_negative::Function)
	
	if any(is_integer.(all_variables(EP)))
		println("The planning model is a MILP")
		optimize!(EP)
			if has_values(EP) #
				planning_output = extractor_planning_problem(EP)
				planning_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in planning_variables]), planning_output = planning_output, theta = value.(EP[:vTHETA])) 
			else
				compute_conflict!(EP)
				list_of_conflicting_constraints = ConstraintRef[];
				for (F, S) in list_of_constraint_types(EP)
					for con in all_constraints(EP, F, S)
						if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
							push!(list_of_conflicting_constraints, con)
						end
					end
				end
				display(list_of_conflicting_constraints)
				@error "The planning solution failed. This should not happen"
			end
	else 
		### The planning model is an LP
		optimize!(EP)
		if has_values(EP)
			neg_cap_bool = check_negative(EP);
			
			if neg_cap_bool
				println("***Resolving the planning problem with Crossover=1 because of negative capacities***")
				set_attribute(EP, "Crossover", 1)
				#set_attribute(EP, "BarHomogeneous", 1)
				optimize!(EP)
				if has_values(EP)
					planning_output = extractor_planning_problem(EP)
					planning_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]), values =Dict([s=>value.(variable_by_name(EP,s)) for s in planning_variables]),planning_output = planning_output, theta = value.(EP[:vTHETA])) 
					set_attribute(EP, "Crossover", 0)
					#set_attribute(EP, "BarHomogeneous", -1)
				else			
					println("The planning problem solution failed, trying with BarHomogenous=1")
					set_attribute(EP, "BarHomogeneous", 1)
					optimize!(EP)
					if has_values(EP)
						planning_output = extractor_planning_problem(EP)
						planning_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]),values =Dict([s=>value.(variable_by_name(EP,s)) for s in planning_variables]), planning_output = planning_output, theta = value.(EP[:vTHETA])) 
						set_attribute(EP, "BarHomogeneous", -1)
					else
						@error "The planning solution failed. This should not happen"
					end
				end
			else
				planning_output = extractor_planning_problem(EP)
				planning_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]),values =Dict([s=>value.(variable_by_name(EP,s)) for s in planning_variables]), planning_output = planning_output, theta = value.(EP[:vTHETA])) 
			end
		else
			println("The planning problem solution failed, trying with BarHomogenous=1")
			set_attribute(EP, "BarHomogeneous", 1)
			optimize!(EP)
			if has_values(EP)
				planning_output = extractor_planning_problem(EP)
				planning_sol =  (LB = objective_value(EP), inv_cost =value(EP[:eObj]),values =Dict([s=>value.(variable_by_name(EP,s)) for s in planning_variables]), planning_output = planning_output, theta = value.(EP[:vTHETA])) 
				set_attribute(EP, "BarHomogeneous", -1)
			else
				@error "The planning solution failed. This should not happen"
			end

		end
	end

	return planning_sol

end

function default_check_negative_capacities(EP::Model)
	neg_cap_bool = false;
	v = all_variables(EP)
	tol = 1e-8
	if minimum(value.(v)) < -tol
		neg_cap_bool = true;
	end
	return neg_cap_bool
end


function solve_dist_subproblems(EP_subproblems::DArray{Dict{Any, Any}, 1, Vector{Dict{Any, Any}}},planning_sol::NamedTuple,inputs, extractor_subproblem::Function)

    p_id = workers();
    np_id = length(p_id);

    sub_results = [Dict() for k in 1:np_id];

    @sync for k in 1:np_id
              @async sub_results[k]= @fetchfrom p_id[k] solve_local_subproblem(localpart(EP_subproblems),planning_sol,inputs, extractor_subproblem); ### This is equivalent to fetch(@spawnat p .....)
    end

	sub_results = merge(sub_results...);

    return sub_results
end

function solve_local_subproblem(subproblem_local::Vector{Dict{Any,Any}},planning_sol::NamedTuple,inputs, extractor_subproblem::Function)

    local_sol=Dict();
    for m in subproblem_local
        EP = m["Model"];
        planning_variables_sub = m["planning_variables_sub"]
        w = m["SubPeriod"];
		local_sol[w] = solve_subproblem(EP,planning_sol,planning_variables_sub,inputs, extractor_subproblem);
    end
    return local_sol
end

function solve_subproblem(EP::Model,planning_sol::NamedTuple,planning_variables_sub::Vector{String},inputs, extractor_subproblem::Function)

	
	fix_planning_variables!(EP,planning_sol,planning_variables_sub)

	optimize!(EP)
	
	if has_values(EP)
		op_cost = objective_value(EP);
        subproblem_results = extractor_subproblem(EP)
		lambda = [dual(FixRef(variable_by_name(EP,y))) for y in planning_variables_sub];
		theta_coeff = 1;
		if haskey(EP,:eObjSlack)
			feasibility_slack = value(EP[:eObjSlack]);
		else
			feasibility_slack = 0.0;
		end
		
	else
		op_cost = 0;
        subproblem_results = Dict{Symbol, Any}();
        compute_conflict!(EP)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(EP)
            for con in all_constraints(EP, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        @info display(list_of_conflicting_constraints)
		lambda = zeros(length(planning_variables_sub));
		theta_coeff = 0;
		feasibility_slack = 0;
		@warn "The subproblem solution failed. This should not happen, double check the input files"
	end
    
	return (op_cost=op_cost,subproblem_results = subproblem_results, lambda = lambda,theta_coeff=theta_coeff,feasibility_slack=feasibility_slack)

end

function fix_planning_variables!(EP::Model,planning_sol::NamedTuple,planning_variables_sub::Vector{String})
	for y in planning_variables_sub
		vy = variable_by_name(EP,y);
		fix(vy,planning_sol.values[y];force=true)
		if is_integer(vy)
			unset_integer(vy)
		elseif is_binary(vy)
			unset_binary(vy)
		end
	end
end
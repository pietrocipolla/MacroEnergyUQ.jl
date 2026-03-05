"""
    run_benders(master_factory, subproblem_factory, distributed, write_outputs, master_optimizer; kwargs...)
    Run Benders decomposition to find optimal solution. 


"""

function run_benders(master_factory::Function, # Function to build master problem ----- must have original objective function labeled as :eObj and theta variables labeled as :vTHETA. Must include optimizer builder.
    subproblem_factory::Function, # Function to build subproblems
    case_path::String;
    workers::Int = 1,
    distributed::Bool = false,
    write_outputs::Bool = false,
    outputs_dir::String = "benders_outputs",
    extractor_subproblem::Function = m -> value.(all_variables(m)), # Function to extract subproblem solutions
	extractor_master::Function = m -> value.(all_variables(m)), # Function to extract master solution
    check_negative_capacities::Function = m -> default_check_negative_capacities(m), # Function to check for negative capacities in the planning problem
    kwargs...
)


    # Build models
    planning_problem, planning_variables, inputs, setup = master_factory(case_path; kwargs...)
    subproblems, planning_variables_sub = subproblem_factory(case_path, planning_variables, setup, inputs; kwargs...)

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
        inputs,
        setup;
        extractor_subproblem, # Function to extract subproblem solutions
        extractor_master, # Function to extract master solution
        check_negative_capacities, # Function to check for negative capacities in the planning problem
        kwargs...
    )

    return benders_result

end

#=============================================================================
    Run Benders with Multiple MGA Iterations


=============================================================================#

function run_benders_mga(
    master_factory::Function, # Function to build master problem ----- must have original objective function labeled as :eObj and theta variables labeled as :vTHETA. Must include optimizer builder.
    subproblem_factory::Function, # Function to build subproblems
    case_path::String; 
    ##### Optional arguments with default (non mga)
    workers::Int = 1,
    distributed::Bool = false,
    write_outputs::Bool = false,
    outputs_dir::String = "benders_outputs",
    extractor_subproblem::Function = m -> value.(all_variables(m)), # Function to extract subproblem solutions
	extractor_master::Function = m -> value.(all_variables(m)), # Function to extract master solution
    check_negative_capacities::Function = m -> default_check_negative_capacities(m), # Function to check for negative capacities in the planning problem
    ##### MGA-specific optional arguments
    n_iterations::Int = 1, # Number of MGA iterations to perform
    vector_factory::Function = default_vector_generator, # Function to generate objective vectors for MGA iterations
    objective_factory::Function = default_objective_factory!,
    kwargs...
)
    println("=" ^ 60)
    println("Phase 1: Running initial Benders optimization")
    println("=" ^ 60)
    
    # Check for distributed computing
    if distributed
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
        inputs,
        setup;
        workers = workers,
        distributed = distributed,
        write_outputs = write_outputs,
        outputs_dir = outputs_dir,
        extractor_master = extractor_master,
        extractor_subproblem = extractor_subproblem,
        check_negative_capacities = check_negative_capacities,
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
                extractor_master = extractor_master,
                extractor_subproblem = extractor_subproblem
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




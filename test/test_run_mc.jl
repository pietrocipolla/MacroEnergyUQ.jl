using Test
using MacroEnergyUQ
using JuMP
using HiGHS
using Distributed

@testset "run_mc tests" begin
    
    # Simple model factory for testing
    function simple_model(optimizer; kwargs...)
        model = Model(optimizer)
        @variable(model, x >= 0)
        @variable(model, y >= 0)
        @objective(model, Min, x + y)
        @constraint(model, x + y >= 1)
        set_silent(model)  # Suppress solver output during tests
        return model
    end
    
    @testset "Basic functionality" begin
        # Test with 5 samples
        data = [0.1 0.2 0.3 0.4 0.5;
                0.6 0.7 0.8 0.9 1.0]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer)
        
        @test length(results.status) == 5
        @test length(results.solve_time) == 5
        @test size(results.outputs) == (5, 2)
        @test all(results.solve_time .>= 0)
    end
    
    @testset "Single sample" begin
        data = reshape([0.5, 0.5], 2, 1)
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer)
        
        @test length(results.status) == 1
        @test size(results.outputs) == (1, 2)
    end

    @testset "Supported model factory return values" begin
        data = [0.1 0.2; 0.6 0.7]

        wrapped_model(optimizer; kwargs...) = (model=simple_model(optimizer; kwargs...),)
        wrapped_results = run_mc(wrapped_model, data, ["x", "y"], HiGHS.Optimizer;
                                 distributed=false)
        @test size(wrapped_results.outputs) == (2, 2)
        @test all(wrapped_results.status .== MOI.OPTIMAL)

        context_model(optimizer; kwargs...) = (
            model=simple_model(optimizer; kwargs...),
            context=(label=42.0,),
        )
        context_extract = (model; ctx) -> [objective_value(model), ctx.index, ctx.label]
        context_results = run_mc(context_model, data, ["x", "y"], HiGHS.Optimizer;
                                 extract=context_extract, distributed=false)
        @test context_results.outputs[:, 2] == [1.0, 2.0]
        @test context_results.outputs[:, 3] == [42.0, 42.0]

        invalid_factory(optimizer; kwargs...) = (unexpected=true,)
        @test_throws ArgumentError MacroEnergyUQ.run_cluster(
            1, invalid_factory, data, ["x", "y"], HiGHS.Optimizer;
            clusters=ones(Int, size(data, 2)),
        )
    end
    
    @testset "Custom extract function" begin
        data = [0.1 0.2 0.3;
                0.6 0.7 0.8]
        
        # Extract only the objective value
        custom_extract = m -> [objective_value(m)]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer; 
                        extract=custom_extract)
        
        @test size(results.outputs) == (3, 1)
        # Check that objective values are positive and finite
        @test all(isfinite.(results.outputs[:, 1]))
        @test all(results.outputs[:, 1] .> 0)
    end
    
    @testset "Multiple clusters" begin
        data = [0.1 0.2 0.3 0.4 0.5 0.6;
                0.6 0.7 0.8 0.9 1.0 1.1]
        
        # Assign samples to 3 clusters
        clusters = [1, 1, 2, 2, 3, 3]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer;
                        clusters=clusters)
        
        @test length(results.status) == 6
        @test size(results.outputs) == (6, 2)
    end
    
    @testset "Non-distributed execution" begin
        data = [0.1 0.2 0.3;
                0.6 0.7 0.8]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer;
                        distributed=false)
        
        @test length(results.status) == 3
        @test size(results.outputs) == (3, 2)
    end
    
    @testset "Input validation" begin
        # Mismatched params and data rows
        data = [0.1 0.2; 0.6 0.7; 0.9 1.0]
        @test_throws ErrorException run_mc(simple_model, data, ["x", "y"], 
                                          HiGHS.Optimizer)
        
        # Mismatched clusters and data columns
        data = [0.1 0.2 0.3; 0.6 0.7 0.8]
        clusters = [1, 1]  # Should be length 3
        @test_throws ErrorException run_mc(simple_model, data, ["x", "y"], 
                                          HiGHS.Optimizer; clusters=clusters)
        
        # Invalid params_type
        data = [0.1 0.2; 0.6 0.7]
        @test_throws ArgumentError run_mc(simple_model, data, ["x", "y"], 
                                         HiGHS.Optimizer; params_type=:invalid)
    end
    
    @testset "Extract function with custom expressions" begin
        data = [0.1 0.2 0.3;
                0.6 0.7 0.8]
        
        # Extract multiple custom values
        custom_extract = m -> [
            value(m[:x]),
            value(m[:y]),
            value(m[:x]) + value(m[:y])
        ]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer; 
                        extract=custom_extract)
        
        @test size(results.outputs) == (3, 3)
        # Third column should be sum of first two
        for i in 1:3
            @test results.outputs[i, 3] ≈ results.outputs[i, 1] + results.outputs[i, 2] atol=1e-6
        end
    end
    
    @testset "Solver status check" begin
        data = [0.1 0.2 0.3;
                0.6 0.7 0.8]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer)
        
        # All solutions should be optimal
        @test all(status -> status == MOI.OPTIMAL, results.status)
    end
    
    @testset "Model with kwargs" begin
        function model_with_constraint(optimizer; min_sum=1.0)
            model = Model(optimizer)
            @variable(model, x >= 0)
            @variable(model, y >= 0)
            @objective(model, Min, x + y)
            @constraint(model, x + y >= min_sum)
            set_silent(model)  # Suppress solver output during tests
            return model
        end
        
        data = [0.1 0.2 0.3;
                0.6 0.7 0.8]
        
        results = run_mc(model_with_constraint, data, ["x", "y"], HiGHS.Optimizer;
                        min_sum=2.0)
        
        @test length(results.status) == 3
        @test all(sum(results.outputs[i, :]) >= 2.0 for i in 1:3)
    end
    
    @testset "Output order consistency with input" begin
        # Create data where we can verify the order matters
        data = [1.0 2.0 3.0 4.0 5.0;
                0.5 0.5 0.5 0.5 0.5]
        
        clusters = [1, 1, 2, 2, 3]
        
        results = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer;
                        clusters=clusters,
                        distributed=false)
        
        # Verify outputs are in the same order as input
        # Since objective coefficients are set from data, 
        # sample i should use coefficients from column i of data
        @test length(results.status) == 5
        @test size(results.outputs) == (5, 2)
        
        # Each result should correspond to its input column
        # The objective is min(data[1,i]*x + data[2,i]*y) s.t. x+y >= 1
        # For our data with data[2,:] all equal to 0.5, 
        # the optimal solution should have x at constraint bound
        for i in 1:5
            @test results.outputs[i, 1] + results.outputs[i, 2] ≈ 1.0 atol=1e-6
        end
    end
    
    @testset "Distributed computing with multiple workers" begin
        # Get initial number of workers
        initial_workers = nworkers()
        
        # Add workers for distributed computing
        if nworkers() == 1
            addprocs(2)
        end
        
        # Load required packages on all workers
        @everywhere using MacroEnergyUQ, JuMP, HiGHS
        
        data = [0.1 0.2 0.3 0.4 0.5 0.6;
                0.6 0.7 0.8 0.9 1.0 1.1]
        
        # Assign samples to different clusters for parallel execution
        clusters = [1, 1, 2, 2, 3, 3]
        
        # Run with distributed=true
        results_dist = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer;
                             clusters=clusters,
                             distributed=true)
        
        # Run without distributed for comparison
        results_seq = run_mc(simple_model, data, ["x", "y"], HiGHS.Optimizer;
                            clusters=clusters,
                            distributed=false)
        
        # Results should be the same
        @test length(results_dist.status) == 6
        @test size(results_dist.outputs) == (6, 2)
        @test all(results_dist.status .== MOI.OPTIMAL)
        
        # Check that outputs match between distributed and sequential
        # This verifies that output order is preserved regardless of execution mode
        @test results_dist.outputs ≈ results_seq.outputs atol=1e-6
        @test results_dist.status == results_seq.status
        
        # Clean up workers if we added them
        if nworkers() > initial_workers
            rmprocs(workers()[end-1:end])
        end
    end
end

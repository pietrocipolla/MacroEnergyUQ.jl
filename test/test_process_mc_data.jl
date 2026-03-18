using Test
using MacroEnergyUQ
using QuasiMonteCarlo
using LinearAlgebra
using RCall

@testset "process_mc_data.jl" begin
    @testset "Basic functionality" begin
        # Test setup
        n_dims = 3      # number of dimensions
        n_points = 100  # number of sample points
        n_threads = 4   # number of parallel threads
        
        # Generate test data using Sobol sequences
        data = QuasiMonteCarlo.sample(n_points, n_dims, SobolSample())
        params = fill(0.5, n_dims)  # reference point in the middle of the hypercube
        
        # Run the function
        processed_data, cluster_assignments, original_indices = process_mc_data(data, n_threads; starting_point = params)
        
        # Basic tests
        @test size(processed_data, 1) == n_dims       # Check dimensions preserved
        @test size(processed_data, 2) == n_points     # Check number of points preserved
        @test length(cluster_assignments) == n_points  # Check cluster assignments length
        @test all(1 .<= cluster_assignments .<= n_threads) # Check cluster assignments range
        
        # Check that all original points are present (just reordered)
        @test all(data[:, original_indices] .== processed_data)
        
        # Test cluster sizes are relatively balanced
        cluster_sizes = [count(==(i), cluster_assignments) for i in 1:n_threads]
        expected_size = n_points ÷ n_threads
        @test all(abs.(cluster_sizes .- expected_size) .<= ceil(Int, n_points/n_threads/2))
    end

    @testset "Edge cases" begin
        # Test with single thread
        data_small = QuasiMonteCarlo.sample(10, 2, SobolSample())
        params_small = [0.5, 0.5]
        processed_data_single, clusters_single = process_mc_data(data_small, 1; starting_point = params_small)
        @test all(clusters_single .== 1)
    end

    @testset "Error handling" begin
        data_error = QuasiMonteCarlo.sample(10, 3, SobolSample())
        params_error = [0.5, 0.5]  # wrong dimension
        
        # Test dimension mismatch
        @test_throws ArgumentError process_mc_data(data_error, 2; starting_point = params_error)
        
        # Test invalid number of threads
        @test_throws ArgumentError process_mc_data(data_error, 0; starting_point = [0.5, 0.5, 0.5])
        @test_throws ArgumentError process_mc_data(data_error, -1; starting_point = [0.5, 0.5, 0.5])
    end
end
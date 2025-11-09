using Test
using MacroEnergyUQ
using QuasiMonteCarlo
using LinearAlgebra

# List of test files
test_files = [
    "test_process_mc_data.jl",
    "test_run_mc.jl"
]

# Run all tests
@testset "MacroEnergyUQ.jl" begin
    for test_file in test_files
        include(test_file)
    end
end
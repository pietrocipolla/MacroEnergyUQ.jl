"""
Tutorial: run MacroEnergyUQ on a monolithic MacroEnergy example case.

The script activates its bundled example environment, which tracks MacroEnergy's
`main` branch, and downloads `multisector_3zone_simpleinputs` from the `main`
branch of MacroEnergyExamples.jl when needed.
"""

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.update("MacroEnergy")

using HiGHS
using JuMP
using MacroEnergy
using MacroEnergyUQ
using Random

# Tutorial settings -----------------------------------------------------------

const EXAMPLE_NAME = "multisector_3zone_simpleinputs"
const EXAMPLES_DIR = get(
    ENV,
    "MACROENERGY_EXAMPLES_DIR",
    joinpath(tempdir(), "MacroEnergyUQExamples"),
)
const CASE_PATH = joinpath(EXAMPLES_DIR, EXAMPLE_NAME)
const N_SAMPLES = parse(Int, get(ENV, "MACROENERGY_UQ_SAMPLES", "20"))
const SEED = 42
const VARIATION = 0.20

# Only expandable technologies whose vNEWUNIT name contains one of these terms
# are varied. Set this to `nothing` to vary every vNEWUNIT coefficient.
const PARAMETER_KEYWORDS = (                                      # CONFIGURE
    "BECCS_Electricity",
    "battery",
    "offshorewind",
    "utilitypv",
    "landbasedwind",
)

const OPTIMIZER_ATTRIBUTES = (
    "solver" => "ipm",
    "run_crossover" => "off",
    "ipm_optimality_tolerance" => 1e-3,
    "output_flag" => false,
)

if !isdir(CASE_PATH)
    mkpath(EXAMPLES_DIR)
    MacroEnergy.download_example(EXAMPLE_NAME, EXAMPLES_DIR; branch = "main")
end
isdir(CASE_PATH) || error("MacroEnergy failed to download $(EXAMPLE_NAME).")


# Model adapter ---------------------------------------------------------------

# Build one model and its extraction context for each Monte Carlo cluster.
function macroenergy_factory(optimizer; case_path, optimizer_attributes)
    case = MacroEnergy.load_case(case_path; lazy_load = true)
    for system in case.systems
        system.settings = merge(
            system.settings,
            (EnableJuMPStringNames = true,),
        )
    end
    configured_optimizer = MacroEnergy.create_optimizer(
        optimizer,
        nothing,
        optimizer_attributes,
    )
    algorithm = MacroEnergy.solution_algorithm(case)
    model = MacroEnergy.generate_model(case, configured_optimizer, algorithm)

    return (model = model, context = (case = case,))
end


# Discover uncertain objective coefficients ----------------------------------

function newunit_parameters(model; keywords = nothing)
    names = String[]
    for variable in all_variables(model)
        variable_name = name(variable)
        is_newunit = startswith(variable_name, "vNEWUNIT_")
        matches = keywords === nothing || any(
            keyword -> occursin(lowercase(keyword), lowercase(variable_name)),
            keywords,
        )
        is_newunit && matches && push!(names, variable_name)
    end
    sort!(names)

    isempty(names) && error(
        "No matching vNEWUNIT variables were found. Check PARAMETER_KEYWORDS " *
        "and ensure EnableJuMPStringNames is true in the case settings.",
    )
    return names
end

reference = macroenergy_factory(
    HiGHS.Optimizer;
    case_path = CASE_PATH,
    optimizer_attributes = OPTIMIZER_ATTRIBUTES,
)
parameter_names = newunit_parameters(
    reference.model;
    keywords = PARAMETER_KEYWORDS,
)
objective = objective_function(reference.model)
base_coefficients = [
    coefficient(objective, variable_by_name(reference.model, parameter_name))
    for parameter_name in parameter_names
]

# Each row is a parameter and each column is one Monte Carlo sample. Here every
# investment coefficient is varied independently by ±VARIATION around its
# nominal value.
rng = MersenneTwister(SEED)
multipliers = 1 .+ VARIATION .* (
    2 .* rand(rng, length(parameter_names), N_SAMPLES) .- 1
)
data = base_coefficients .* multipliers

# Reorder nearby samples; one cluster avoids the optional RCall extension.
processed_data, clusters, original_indices = process_mc_data(
    data,
    1;
    starting_point = base_coefficients,
)


# Solve and collect results ---------------------------------------------------

# Extract the objective and selected capacity variables from each solution.
function extract_summary(model; ctx)
    capacities = [
        value(variable_by_name(model, parameter_name))
        for parameter_name in parameter_names
    ]
    return [objective_value(model); capacities]
end

# Update objective coefficients and solve every sample.
results = run_mc(
    macroenergy_factory,
    processed_data,
    parameter_names,
    HiGHS.Optimizer;
    params_type = :objective,
    clusters,
    extract = extract_summary,
    optimizer_attributes = OPTIMIZER_ATTRIBUTES,
    case_path = CASE_PATH,
)

# `process_mc_data` changes sample order; restore the original generated order.
outputs = similar(results.outputs)
statuses = similar(results.status)
solve_times = similar(results.solve_time)
outputs[original_indices, :] = results.outputs
statuses[original_indices] = results.status
solve_times[original_indices] = results.solve_time

output_names = ["objective"; parameter_names]
@info "MacroEnergy uncertainty run complete" N_SAMPLES parameter_names output_names statuses outputs

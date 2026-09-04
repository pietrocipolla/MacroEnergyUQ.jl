"""
Template: run MacroEnergyUQ on a monolithic MacroEnergy case.

Adapt every value marked CONFIGURE. Parameter names are scalar JuMP variable
names whose objective coefficients will be replaced for each sample.
"""

using HiGHS
using JuMP
using MacroEnergy
using MacroEnergyUQ

const CASE_PATH = "/absolute/path/to/macroenergy_case"             # CONFIGURE
const PARAMETER_NAMES = [                                         # CONFIGURE
    "vNEWUNIT_example_1_period1",
    "vNEWUNIT_example_2_period1",
]
const OPTIMIZER_ATTRIBUTES = ("output_flag" => false,)
const WRITE_MACROENERGY_OUTPUTS = false
const N_CLUSTERS = 1

# Rows match PARAMETER_NAMES; columns are samples.                  # CONFIGURE
data = [100.0 110.0 120.0;
        150.0 145.0 140.0]

function macroenergy_factory(
    optimizer;
    case_path,
    optimizer_attributes,
)
    case = MacroEnergy.load_case(case_path; lazy_load = true)
    configured_optimizer = MacroEnergy.create_optimizer(
        optimizer,
        nothing,
        optimizer_attributes,
    )
    model = MacroEnergy.generate_model(case, configured_optimizer)

    # MacroEnergyUQ adds the sample index to this context before extraction.
    return (model, context = (case = case,))
end


function macroenergy_extract(model; ctx)
    # CONFIGURE: append scalar JuMP values needed by the analysis, e.g.
    # push!(output, value(variable_by_name(model, "vNEWUNIT_..._period1")))
    output = [objective_value(model)]

    if WRITE_MACROENERGY_OUTPUTS
        sample_dir = joinpath(
            @__DIR__,
            "macroenergy_outputs",
            "sample_$(ctx.index)",
        )
        mkpath(sample_dir)
        MacroEnergy.write_outputs(sample_dir, ctx.case, model)
    end

    return output
end


# Start with one cluster, which does not require the optional RCall extension.
# For multiple clusters, load RCall and use process_mc_data with N_CLUSTERS.
N_CLUSTERS == 1 || error(
    "Load RCall and adapt this template before using multiple clusters.",
)
processed_data, clusters, original_indices = process_mc_data(data, N_CLUSTERS)

results = run_mc(
    macroenergy_factory,
    processed_data,
    PARAMETER_NAMES,
    HiGHS.Optimizer;
    params_type = :objective,
    clusters,
    extract = macroenergy_extract,
    optimizer_attributes = OPTIMIZER_ATTRIBUTES,
    case_path = CASE_PATH,
)

# Map the processed result rows back to the original sample order.
outputs_in_original_order = similar(results.outputs)
outputs_in_original_order[original_indices, :] = results.outputs

@info "MacroEnergy uncertainty run complete" results.status outputs_in_original_order

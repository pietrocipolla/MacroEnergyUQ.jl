"""
Template: run MacroEnergyUQ on a single-stage GenX case.

Adapt every value marked CONFIGURE. Parameter names are scalar JuMP variable
names whose objective coefficients will be replaced for each sample.
"""

using GenX
using HiGHS
using JuMP
using MacroEnergyUQ

const CASE_PATH = "/absolute/path/to/genx_case"                  # CONFIGURE
const PARAMETER_NAMES = ["vCAP[1]", "vCAP[2]"]                 # CONFIGURE
const WRITE_GENX_OUTPUTS = false
const N_CLUSTERS = 1

# Rows match PARAMETER_NAMES; columns are samples.                  # CONFIGURE
data = [100.0 110.0 120.0;
        150.0 145.0 140.0]

function genx_factory(optimizer; case_path)
    settings_path = GenX.get_settings_path(case_path)
    setup = GenX.configure_settings(
        joinpath(settings_path, "genx_settings.yml"),
        joinpath(settings_path, "output_settings.yml"),
    )
    setup["MultiStage"] == 0 || error(
        "This template supports single-stage GenX cases only.",
    )

    configured_optimizer = GenX.configure_solver(settings_path, optimizer)
    inputs = GenX.load_inputs(setup, case_path)
    model = GenX.generate_model(setup, inputs, configured_optimizer)

    # Context is available to genx_extract as `ctx`; MacroEnergyUQ adds ctx.index.
    return (model, context = (setup = setup, inputs = inputs))
end

function genx_extract(model; ctx)
    # CONFIGURE: add scalar outputs required by the analysis. JuMP containers
    # can be indexed directly, for example model[:vCAP][resource_id].
    output = [objective_value(model)]

    if WRITE_GENX_OUTPUTS
        sample_dir = joinpath(@__DIR__, "genx_outputs", "sample_$(ctx.index)")
        mkpath(sample_dir)
        ctx.inputs["solve_time"] = solve_time(model)
        GenX.write_outputs(model, sample_dir, ctx.setup, ctx.inputs)
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
    genx_factory,
    processed_data,
    PARAMETER_NAMES,
    HiGHS.Optimizer;
    params_type = :objective,
    clusters,
    extract = genx_extract,
    case_path = CASE_PATH,
)

# Map the processed result rows back to the original sample order.
outputs_in_original_order = similar(results.outputs)
outputs_in_original_order[original_indices, :] = results.outputs

@info "GenX uncertainty run complete" results.status outputs_in_original_order

# [MacroEnergy.jl tutorial](@id macroenergy-tutorial)

This tutorial runs a monolithic MacroEnergy case over uncertain investment
costs. It discovers the relevant capacity variables, samples their objective
coefficients, solves every sample, and collects objectives and capacities.

The model uses the `main` branch of
[`MacroEnergy.jl`](https://github.com/macroenergy/MacroEnergy.jl). The case is
the `multisector_3zone_simpleinputs` example created and maintained by the
MacroEnergy team in
[`MacroEnergyExamples.jl`](https://github.com/macroenergy/MacroEnergyExamples.jl).
The case is downloaded from that repository's `main` branch; MacroEnergyUQ does
not copy or modify the upstream data.

The complete script is available as
[`examples/macroenergy_template.jl`](https://github.com/pietrocipolla/MacroEnergyUQ.jl/blob/master/examples/macroenergy_template.jl).
The sections below unpack that script one step at a time.

## 1. Prepare the environment and example case

The standalone script activates the included `examples/Project.toml` and
instantiates it, so no separate project has to be created. That project obtains
MacroEnergy directly from its GitHub `main` branch. On its first run the script
also downloads `multisector_3zone_simpleinputs` from the examples repository's
`main` branch through MacroEnergy's `download_example` function. The data are
cached under the system temporary directory and reused on later runs.

Set `MACROENERGY_EXAMPLES_DIR` to choose a different cache location. No case
path needs to be supplied.

`PARAMETER_KEYWORDS` limits the study to technologies whose `vNEWUNIT_...`
variable name contains one of the listed terms. Set it to `nothing` to include
every `vNEWUNIT_...` variable.

```julia
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.update("MacroEnergy")

using HiGHS
using JuMP
using MacroEnergy
using MacroEnergyUQ
using Random

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

const PARAMETER_KEYWORDS = (
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
```

## 2. Define the model factory

[`run_mc`](@ref) receives a *model factory*: a function that constructs a clean
model when given an optimizer. Extra keyword arguments passed to `run_mc` are
forwarded to this function.

MacroEnergyUQ calls the factory once per cluster and reuses the resulting model
for every sample in that cluster. The factory below loads the case, enables
JuMP string names in memory, configures HiGHS through MacroEnergy, selects the
case's configured solution algorithm, and generates the JuMP model. Enabling
names in memory leaves the downloaded example unchanged.

The return value contains both `model` and `context`. Context is optional, but
it is useful when extraction needs MacroEnergy-specific objects. Before each
extraction, MacroEnergyUQ adds the current processed-sample index as
`ctx.index`.

```julia
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
```

## 3. Identify the uncertain coefficients

In this example, uncertainty applies to investment costs. A cost is represented
by the objective coefficient of a scalar `vNEWUNIT_...` capacity variable.
`newunit_parameters` selects those variables by name and optionally filters
them by technology keyword.

One reference model is then built to obtain the exact variable names and their
nominal objective coefficients. These names must also exist in every model
subsequently created by the factory.

```julia
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
```

## 4. Generate Monte Carlo samples

MacroEnergyUQ expects parameters along the rows of `data` and samples along its
columns. The code below draws an independent multiplier for every
parameter-sample pair. With `VARIATION = 0.20`, each nominal coefficient is
multiplied by a value uniformly distributed between 0.8 and 1.2.

`MersenneTwister(SEED)` makes the generated study reproducible without changing
Julia's global random-number generator.

```julia
rng = MersenneTwister(SEED)
multipliers = 1 .+ VARIATION .* (
    2 .* rand(rng, length(parameter_names), N_SAMPLES) .- 1
)
data = base_coefficients .* multipliers
```

## 5. Arrange samples for solution

[`process_mc_data`](@ref) places nearby samples next to one another. Because a
cluster reuses the same model, this ordering can help consecutive solves.

It returns the reordered data, one cluster assignment per reordered sample,
and `original_indices`, which records where each reordered sample appeared in
the original data. Using one cluster keeps this tutorial sequential and avoids
the optional RCall extension required for multi-cluster preprocessing.

The nominal coefficients are passed as `starting_point`, so the ordering starts
near the unperturbed case.

```julia
processed_data, clusters, original_indices = process_mc_data(
    data,
    1;
    starting_point = base_coefficients,
)
```

## 6. Define the extraction function

MacroEnergyUQ calls `extract_summary` after each model solve. The function must
return a numeric vector of the same length for every sample; those vectors
become the rows of `results.outputs`.

Here, the first output is the optimized objective and the remaining outputs are
the installed capacities associated with the uncertain investment costs. Since
the factory returned context, the extractor must accept the `ctx` keyword. This
example does not need it, but `ctx.case` and `ctx.index` are available for
case-specific output processing.

```julia
function extract_summary(model; ctx)
    capacities = [
        value(variable_by_name(model, parameter_name))
        for parameter_name in parameter_names
    ]
    return [objective_value(model); capacities]
end
```

## 7. Run the study

For every column of `processed_data`, [`run_mc`](@ref) replaces the objective
coefficient of each variable in `parameter_names`, optimizes the model, invokes
`extract_summary`, and records the termination status and solve time.

`params_type = :objective` selects JuMP's objective-coefficient update path.
The optimizer attributes and case path are not interpreted by MacroEnergyUQ;
they are forwarded to `macroenergy_factory`.

```julia
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
```

## 8. Restore the generated sample order

The rows returned by `run_mc` follow the reordered data. Applying
`original_indices` restores outputs, statuses, and solve times to the order in
which the samples were originally generated.

The first column of `outputs` is the objective. The remaining columns follow
`parameter_names`, as recorded in `output_names`.

```julia
outputs = similar(results.outputs)
statuses = similar(results.status)
solve_times = similar(results.solve_time)
outputs[original_indices, :] = results.outputs
statuses[original_indices] = results.status
solve_times[original_indices] = results.solve_time

output_names = ["objective"; parameter_names]
@info "MacroEnergy uncertainty run complete" N_SAMPLES parameter_names output_names statuses outputs
```

## Run the complete file

Run the script directly from the MacroEnergyUQ repository. It activates and
instantiates the bundled example environment itself:

```bash
julia examples/macroenergy_template.jl
```

Set `MACROENERGY_UQ_SAMPLES` to override the default of 20 samples without
editing the file. Internet access is required to check MacroEnergy's configured
GitHub `main` branch. Installed package files and the downloaded case are
cached locally.

For a Benders MacroEnergy run, use the component return format described in
[Benders factories](@ref) instead of the monolithic `(model, context)` return.

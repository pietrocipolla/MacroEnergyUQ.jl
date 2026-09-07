# Monte Carlo simulations

## Inputs and parameter updates

[`run_mc`](@ref) accepts `data::Matrix{Float64}` with shape
`(number_of_parameters, number_of_samples)`. `params[i]` names the JuMP object
updated from `data[i, sample]`.

Two update modes are supported:

- `params_type = :objective` changes the objective coefficient of the named
  variable with `JuMP.set_objective_coefficient`.
- `params_type = :parameter` changes the value of a named JuMP parameter with
  `JuMP.set_parameter_value`.

!!! important "String names"
    Names must be Julia/JuMP string names, including indices where applicable, for
    example `"vCAP[3]"`. String names must be enabled in `JuMP`. We are planning
    to extend integration for parameters without string names in the future.

## Model factories

A factory is called once per cluster as `factory(optimizer; kwargs...)`. It may
return:

1. a `JuMP.Model`;
2. `(model = model,)`, optionally with `context = (...)`; or
3. Benders components with `planning_problem`, `subproblems`, and
   `linking_variables`, optionally with `settings` and `context`. The code underlying 
   uses [`MacroEnergySolvers.jl`](https://github.com/macroenergy/MacroEnergySolvers.jl).
   Therefore, the components should be compliant with the required structure.

Context must be a `NamedTuple`. Before extraction, MacroEnergyUQ adds the
one-based sample index and calls `extract(model_or_results; ctx)`. This is useful
when an upstream package needs its input data to write native results:

```julia
function factory(optimizer; case_path)
    inputs = load_case_inputs(case_path)
    model = build_model(inputs, optimizer)
    return (model, context = (inputs = inputs, case_path = case_path))
end

function extract(model; ctx)
    sample_dir = joinpath("outputs", "sample_$(ctx.index)")
    write_native_outputs(sample_dir, model, ctx.inputs)
    return [objective_value(model)]
end
```

The extraction function must return a numeric vector of the same length for
every sample. When no extractor is supplied, all JuMP variable values are
returned in `all_variables(model)` order.

## Clusters and parallel execution

`clusters` assigns each sample column to a positive integer cluster. Samples in
one cluster are solved sequentially on the same model; clusters may run in
parallel. For predictable execution, use contiguous labels `1:n` and ensure
cluster `1` is nonempty.

With `distributed = false` (the default), `pmap(...; distributed=false)` uses
Julia task-based execution. With `distributed = true`, add workers and load all
packages and top-level factory/extractor definitions on them before calling
`run_mc`:

```julia
using Distributed
addprocs(4; exeflags = "--project=.")
@everywhere using MacroEnergyUQ, JuMP, HiGHS

# Define factory and extractor with @everywhere at top level.
results = run_mc(
    factory,
    data,
    params,
    HiGHS.Optimizer;
    clusters,
    distributed = true,
)
```

Use `cluster_start_delay_s` to stagger cluster startup when model creation
contends for shared files or licences. Results are always restored to input
sample order.

## Writing outputs to disk

For a small fixed-length output vector, keep the default `write_outputs=false`.
For large vectors:

```julia
results = run_mc(
    factory,
    data,
    params,
    optimizer;
    write_outputs = true,
    output_dir = "mc_outputs",
)
```

Each vector is written as `sample_<index>.csv`; the return value contains
`status`, `solve_time`, and `output_dir` rather than `outputs`. If an upstream
package writes a directory of native outputs, do that inside an extractor with
context and keep `write_outputs=false` for a compact summary vector.

## Benders factories

For a MacroEnergySolvers workflow, return:

```julia
(
    planning_problem = planning_problem,
    subproblems = subproblems,
    linking_variables = linking_variables,
    settings = benders_settings, # optional; defaults to Dict()
    context = (case = case,),    # optional
)
```

Uncertain parameters are updated on `planning_problem`. MacroEnergyUQ calls
`MacroEnergySolvers.benders` for every sample and passes its result to the
extractor. Benders statuses are currently reported as the string `"BENDERS"`.

## Failure handling

`run_mc` does not catch solver or extractor exceptions. Check `results.status`
before interpreting outputs, and make custom extractors robust to the statuses
your model can produce. A failed cluster aborts the call rather than returning
a partial result.

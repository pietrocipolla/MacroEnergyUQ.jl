# Getting started

## Installation

Install the package directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/pietrocipolla/MacroEnergyUQ.jl")
```

For development from a local checkout, use `Pkg.develop(path = "...")`.

## A complete small example

The factory receives the optimizer as its first positional argument. It must
construct a fresh model because each cluster is processed independently.

```julia
using HiGHS
using JuMP
using MacroEnergyUQ

function capacity_model(optimizer; demand = 100.0)
    model = Model(optimizer)
    set_silent(model)

    @variable(model, wind >= 0)
    @variable(model, gas >= 0)
    @constraint(model, wind + gas >= demand)
    @objective(model, Min, 40.0 * wind + 70.0 * gas)
    return model
end

# Rows: wind and gas objective coefficients.
# Columns: four Monte Carlo samples.
data = [35.0 45.0 55.0 65.0;
        75.0 70.0 65.0 60.0]

extract(model) = [
    objective_value(model),
    value(model[:wind]),
    value(model[:gas]),
]

results = run_mc(
    capacity_model,
    data,
    ["wind", "gas"],
    HiGHS.Optimizer;
    extract,
)
```

`results.status` and `results.solve_time` follow the original sample order.
When results are kept in memory, `results.outputs[i, :]` is the vector returned
by `extract` for sample `i`.

## Next steps

- Read [Monte Carlo simulations](@ref) for model factory return types,
  extraction contexts, parallel modes, and disk output.
- Read [Sample preprocessing](@ref) before generating cluster assignments.
- Use the [GenX.jl template](@ref) or [MacroEnergy.jl template](@ref) as a
  case-specific starting point.

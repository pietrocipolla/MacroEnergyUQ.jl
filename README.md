# MacroEnergyUQ.jl

[![codecov](https://codecov.io/github/pietrocipolla/MacroEnergyUQ.jl/branch/master/graph/badge.svg)](https://app.codecov.io/github/pietrocipolla/MacroEnergyUQ.jl)

MacroEnergyUQ is a Julia package for uncertainty quantification of JuMP-based
optimization models, with a focus on energy-system planning. It provides Monte
Carlo execution, sample clustering and reordering, disk-backed outputs, Benders
integration through MacroEnergySolvers, and optional optimal-transport
sensitivity analysis through RCall.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/pietrocipolla/MacroEnergyUQ.jl")
```

## Quick start

```julia
using HiGHS, JuMP, MacroEnergyUQ

function create_model(optimizer; demand = 1.0)
    model = Model(optimizer)
    @variable(model, wind >= 0)
    @variable(model, gas >= 0)
    @constraint(model, wind + gas >= demand)
    @objective(model, Min, 40.0 * wind + 70.0 * gas)
    return model
end

# Each row is a parameter; each column is a sample.
data = [35.0 45.0 55.0;
        75.0 70.0 65.0]

results = run_mc(
    create_model,
    data,
    ["wind", "gas"],
    HiGHS.Optimizer;
    extract = model -> [
        objective_value(model),
        value(model[:wind]),
        value(model[:gas]),
    ],
)
```

`results.outputs` has one row per sample. `results.status` and
`results.solve_time` use the same order as the input columns.

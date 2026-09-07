# MacroEnergyUQ.jl

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

## Documentation and templates

The complete manual is under [`docs/src`](docs/src/index.md), including:

- [getting started](docs/src/getting_started.md);
- [Monte Carlo factories, parallelism, and outputs](docs/src/monte_carlo.md);
- [sample preprocessing](docs/src/preprocessing.md);
- [optimal-transport sensitivity analysis](docs/src/sensitivity.md); and
- the [API reference](docs/src/api.md).

Two case-oriented templates are included:

- [`examples/genx_template.jl`](examples/genx_template.jl)
- [`examples/macroenergy_template.jl`](examples/macroenergy_template.jl)

Both templates mark the case-specific paths, JuMP variable names, sample data,
and extracted outputs that must be adapted.

## Building the manual

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Deploying the manual

The documentation workflow builds and deploys the manual to the `gh-pages`
branch after a push to `master`, a `v*` tag, or a manual run from the GitHub
Actions page. In the repository's **Settings → Pages** screen, select
**Deploy from a branch**, then choose the `gh-pages` branch and `/ (root)`.

The published development documentation is available at
<https://pietrocipolla.github.io/MacroEnergyUQ.jl/dev/>.

RCall is optional for the basic single-cluster Monte Carlo workflow. Load it
when using multi-cluster preprocessing, multivariate quantile transformation,
or the sensitivity-analysis functions.

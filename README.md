# MacroEnergyUQ.jl

A Julia package for uncertainty quantification of optimization models, with a specific focus on energy planning models.

Requires Julia 1.10 or later.

## Features

- Data preprocessing with multivariate quantile transformation
- Integrated clustering for optimized parallel execution
- Automatic parallelization of Monte Carlo simulations
- Support for JuMP models with and without Parameters

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/pietrocipolla/MacroEnergyUQ.jl")
```

### Troubleshooting RCall Installation

If you encounter issues with the RCall installation (which is required for this package), you can try the CondaPkg approach described in the [RCall installation guide](https://juliainterop.github.io/RCall.jl/stable/installation/#(Experimental)-Usage-with-CondaPkg). This package uses CondaPkg to manage the R environment automatically, but in some cases manual configuration may be needed.

## Basic Usage

```julia
using MacroEnergyUQ
using JuMP
using HiGHS
using QuasiMonteCarlo

# Define a function that creates your model
function create_model(optimizer; demand::Float64 = 1.0, min_production::Float64 = 0.0)
    model = Model(optimizer)
    @variable(model, x >= min_production)
    @variable(model, y >= min_production)
    @objective(model, Min, 0.1*x + 0.6*y)
    @constraint(model, x + y >= demand)
    return model
end

# Generate sample points using QuasiMonteCarlo
n_samples = 100
n_params = 2
data = QuasiMonteCarlo.sample(n_samples, n_params, SobolSample())

# Run Monte Carlo simulations
results = run_mc(create_model, data, ["x", "y"], HiGHS.Optimizer;
                demand = 2.0, min_production = 0.1)
```

The model factory must accept the optimizer as its first argument. It may return
the `JuMP.Model` directly, as above, or return `(model=model, context=context)`
when the extraction function needs additional context. In the latter case,
`context` must be a NamedTuple and the extractor is called as
`extract(model; ctx=context)`, with the current one-based sample `index` added
to `context` automatically.

## Advanced Features

### Data Preprocessing

```julia
using Distributions

data = cat(rand(Normal(100, 10), 100)', rand(LogNormal(4, 0.3), 100)', dims=1) 
data, clusters, original_indices = MacroEnergyUQ.process_mc_data(data, 4)
```

### Custom Parallelization

```julia
# Use clustering for optimized parallel execution
results = run_mc(create_model, processed_data, ["x", "y"], Gurobi.Optimizer;
                clusters = clusters,  # Custom cluster assignment
                demand = 2.0,
                min_production = 0.1)
```

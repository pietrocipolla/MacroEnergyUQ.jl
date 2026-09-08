# Sample preprocessing

[`process_mc_data`](@ref) groups samples into clusters and reorders the samples
inside each cluster so consecutive model solves are close in parameter space.

```julia
processed, clusters, original_indices = process_mc_data(
    data,
    4;
    starting_point = fill(0.5, size(data, 1)),
)

results = run_mc(factory, processed, params, optimizer; clusters)
```

As with [`run_mc`](@ref), columns are samples. `processed` is exactly
`data[:, original_indices]`. Consequently, row `i` of `results.outputs`
corresponds to column `original_indices[i]` of the original data. To restore
outputs to original order:

```julia
outputs_in_original_order = similar(results.outputs)
outputs_in_original_order[original_indices, :] = results.outputs
```

## Options

- `n_threads` determines the number of returned clusters.
- `starting_point` is the reference vector for nearest-neighbour ordering. Its
  length must equal the number of parameters.
- `sorting_algorithm = :nearest_neighbor` uses a fast greedy solution to the 
  shortest path problem by using nearest neighbors.
- `sorting_algorithm = :tsp` solves a path optimization problem and requires a
  supplied JuMP `optimizer`. It is the actual solution to the shortest path problem
  but it may be computationally expensive for large datasets.
- `quantile_transform = true` first maps samples to a multivariate uniform
  quantile space using optimal transport. It avoids the problem of having very
  different, possibly correlated, parameter spaces.

The single-cluster nearest-neighbour path works without RCall. Creating more
than one cluster or enabling the quantile transform requires loading RCall so
the optional extension is active:

```julia
using RCall
using MacroEnergyUQ
```

The R environment must provide the packages used by the extension. The package
attempts to configure them through CondaPkg, but system-specific RCall setup may
still be necessary.

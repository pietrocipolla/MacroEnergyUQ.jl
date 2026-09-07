# Sensitivity analysis

MacroEnergyUQ exposes wrappers around the R package `gsaot`. Load RCall before
calling these functions so the package extension is activated:

```julia
using MacroEnergyUQ
using RCall
```

The Monte Carlo convention is `parameters × samples`, while the sensitivity
wrappers transpose `x` internally before sending it to R. Therefore, pass the
same parameter-by-sample input matrix used by [`run_mc`](@ref). The output
matrix returned by `run_mc` already has shape `samples × outputs`:

```julia
indices = ot_indices(data, results.outputs, 10)
threshold = irrelevance_threshold(results.outputs, 10)
```

## Optimal-transport indices

```julia
ot_indices(x, y, M;
    cost = "L2",
    discrete_out = false,
    solver = "sinkhorn",
    solver_optns = nothing,
    scaling = true,
    boot = false,
    stratified_boot = true,
    R = nothing,
    parallel = "no",
    ncpus = 1,
    conf = 0.95,
    type = "norm",
)
```

`M` is the number of partitions for continuous inputs. Enable `boot` and set
`R` to request bootstrap replicas and confidence intervals. The exact returned
dictionary is produced by the installed `gsaot` version.

## Wasserstein--Bures indices

Use `ot_indices_wb(x, y, M; ...)` for the Wasserstein--Bures variant. Its common
bootstrap options are `boot`, `R`, `parallel`, `ncpus`, `conf`, and `type`.

## Irrelevance threshold

`irrelevance_threshold(y, M; ...)` estimates a dummy-variable threshold against
which sensitivity indices can be compared. It accepts `dummy_optns`, `cost`,
`discrete_out`, `solver`, `solver_optns`, and `scaling`.

!!! note
    These functions delegate computation and validation to R. Pin and record
    your R and `gsaot` versions for reproducible studies.

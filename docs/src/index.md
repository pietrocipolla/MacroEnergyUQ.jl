# MacroEnergyUQ.jl

MacroEnergyUQ is a Julia package for uncertainty quantification of JuMP-based
optimization models, with utilities aimed at energy-system planning workflows
with [GenX.jl](https://github.com/GenXProject/GenX.jl) and 
[MacroEnergy.jl](https://github.com/macroenergy/MacroEnergy.jl).
MacroEnergyUQ can run a model repeatedly over a matrix of uncertain inputs, 
and it is specifically designed for large-scale models. It can:
- cluster the inputs for multithreaded and distributed computing using 
  [Distributed.jl](https://docs.julialang.org/en/v1/manual/parallel-computing/);
- reuse one model per sample cluster to reduce model generation overhead;
- reorder and cluster samples to reduce the distance between successive solves, 
  aiding warmstartable algorithms;
- write large per-sample outputs to disk; and
- calculate optimal-transport sensitivity indices through an optional RCall
  extension.

The central function is [`run_mc`](@ref). A user supplies a model factory, a
parameter-by-sample matrix, the matching JuMP variable names, and an optimizer.

```julia
results = run_mc(factory, data, parameter_names, optimizer)
```

!!! important "Data orientation"
    `data` has one uncertain parameter per **row** and one Monte Carlo sample per
    **column**. The output matrix has one sample per row.

Start with [Getting started](@ref), then adapt either the
[GenX.jl template](@ref) or
[MacroEnergy.jl tutorial](@ref macroenergy-tutorial) for a full energy model.

## Package requirements

- Julia 1.10 or later.
- A JuMP-compatible optimizer.
- RCall and the required R packages only when using multi-cluster preprocessing,
  multivariate quantile transformation, or sensitivity analysis.

## Contents

```@contents
Pages = [
    "getting_started.md",
    "monte_carlo.md",
    "preprocessing.md",
    "sensitivity.md",
    "examples/genx.md",
    "examples/macroenergy.md",
    "api.md",
]
Depth = 2
```

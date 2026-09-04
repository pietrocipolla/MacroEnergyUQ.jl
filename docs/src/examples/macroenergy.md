# MacroEnergy.jl template

The repository includes
[`examples/macroenergy_template.jl`](https://github.com/pietrocipolla/MacroEnergyUQ.jl/blob/master/examples/macroenergy_template.jl),
a monolithic MacroEnergy starting point.

It loads a case lazily, creates the optimizer through MacroEnergy, builds a new
model for each cluster, and returns the case as extraction context. The
extractor shows how to collect a compact summary and optionally call
`MacroEnergy.write_outputs` in a sample-specific directory.

Replace `CASE_PATH`, `PARAMETER_NAMES`, `OPTIMIZER_ATTRIBUTES`, and the sample
matrix. The parameter names must match scalar JuMP variable names in the
generated planning model, commonly case-specific `vNEWUNIT_..._period...`
names when varying investment objective coefficients.

```bash
julia --project=. examples/macroenergy_template.jl
```

For a Benders MacroEnergy run, use the component return format described in
[Benders factories](@ref) instead of the monolithic `(model, context)` return.

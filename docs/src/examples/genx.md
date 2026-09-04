# GenX.jl template

The repository includes [`examples/genx_template.jl`](https://github.com/pietrocipolla/MacroEnergyUQ.jl/blob/master/examples/genx_template.jl),
a single-stage GenX starting point.

The template reads case settings, loads GenX inputs, builds a fresh JuMP model
per cluster, and returns inputs/settings as extraction context. It can also
write native GenX output under one directory per sample.

Before running it, replace `CASE_PATH`, `PARAMETER_NAMES`, the sample matrix,
and the extraction variables. A GenX parameter name must be the exact JuMP
string name of a scalar variable such as `"vCAP[3]"`, not merely a CSV column
name. GenX must have `EnableJuMPStringNames` enabled.

```bash
julia --project=. examples/genx_template.jl
```

The template deliberately rejects multi-stage cases because a multi-stage
factory requires a study-specific decision about which stage models and
coefficients are uncertain.

using Documenter

# Load the package checkout rather than a registered release.
pushfirst!(LOAD_PATH, dirname(@__DIR__))
using MacroEnergyUQ

makedocs(
    sitename = "MacroEnergyUQ.jl",
    modules = [MacroEnergyUQ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://pietrocipolla.github.io/MacroEnergyUQ.jl",
        edit_link = nothing,
        repolink = "https://github.com/pietrocipolla/MacroEnergyUQ.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "User guide" => [
            "Monte Carlo simulations" => "monte_carlo.md",
            "Sample preprocessing" => "preprocessing.md",
            "Sensitivity analysis" => "sensitivity.md",
        ],
        "Examples" => [
            "GenX.jl template" => "examples/genx.md",
            "MacroEnergy.jl tutorial" => "examples/macroenergy.md",
        ],
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
    remotes = nothing,
)

deploydocs(
    repo = "github.com/pietrocipolla/MacroEnergyUQ.jl.git",
    devbranch = "master",
)

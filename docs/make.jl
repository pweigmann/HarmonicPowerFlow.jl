using Documenter 
using HarmonicPowerFlow

makedocs(
    modules = [HarmonicPowerFlow],
    authors = "Pascal Weigmann",
    sitename ="HarmonicPowerFlow.jl",
    pages = [
        "General" => "index.md"
        "Conventions" => "conventions.md"
    ])

deploydocs(repo = "github.com/pweigmann/HarmonicPowerFlow.jl.git")
using Documenter
using HarmonicPowerFlow

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [HarmonicPowerFlow],
    authors = "Pascal Weigmann",
    sitename ="HarmonicPowerFlow.jl",
    pages = [
        "General" => "index.md"
        "Conventions" => "conventions.md"
        "Functions" => "functions.md"
    ])

deploydocs(
    repo = "github.com/pweigmann/HarmonicPowerFlow.jl.git",
    devbranch = "main"
    )
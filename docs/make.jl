push!(LOAD_PATH, ".../src/")

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

deploydocs(
    repo = "github.com/pweigmann/HarmonicPowerFlow.jl.git",
    deploy_config = Documenter.GitHubActions("github.com/pweigmann/HarmonicPowerFlow.jl.git", "push", "refs/heads/master")
)
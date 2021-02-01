include("../src/MCPhylo.jl")
using .MCPhylo
using Documenter

makedocs(root="./",
    source  = "src",
    build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[MCPhylo],
    sitename="MCPhylo",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "intro.md",
        "Tree Functionality" => "Tree.md",
        "Likelihood" => "Likelihood.md",
        "Parser" => "Parser.md",
        "Distributions" => "distributions.md",
        "Model" => "model.md",
        "Output" => "output.md",
        "Samplers" => "samplers.md",
        "Utils" => "Utils.md",
        "Links" => "Links.md"]
   )

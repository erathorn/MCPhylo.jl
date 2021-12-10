using Documenter
using MCPhylo

makedocs(
    modules = MCPhylo,
    sitename="MCPhylo.jl",
    format = Documenter.HTML(),
    pages = [
        "Index" => "index.md",
        "Likelihood" => "likelihood.md",
        "Distributions" => "distributions.md",
        "Model" => "model.md",
        "Output" => "output.md",
        "Samplers" => "samplers.md",
        "Utils" => "utils.md",
        "Links" => "links.md"]
   )

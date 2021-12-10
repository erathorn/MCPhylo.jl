using Documenter
using MCPhylo

makedocs(
    modules = MCPhylo,
    sitename="MCPhylo.jl",
    format = Documenter.HTML(),
    pages = [
        "Index" => "index.md",
        "Likelihood" => "likelihood1.md",
        "Distributions" => "distributions.md",
        "Model" => "model.md",
        "Output" => "output.md",
        "Samplers" => "samplers.md",
        "Utils" => "utils.md",
        "Links" => "links.md"]
   )

deploydocs(
    repo = "github.com/erathorn/MCPhylo.jl.git",
    devbranch="main"
     )


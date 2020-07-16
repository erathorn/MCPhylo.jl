include("../src/MCPhylo.jl")
using .MCPhylo
using Documenter

makedocs(sitename="MCPhylo",
    modules = [MCPhylo],
    format = Documenter.HTML(prettyurls = false)
   )

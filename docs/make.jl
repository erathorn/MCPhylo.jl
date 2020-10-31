include("../src/MCPhylo.jl")
using .MCPhylo
using Documenter

makedocs(root="./",
    source  = "src",
    build   = "build",
    clean   = true,
    doctest = true,
    modules = [MCPhylo],
    sitename="MCPhylo",
    format = Documenter.HTML(prettyurls = false)
   )

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
        "Adhams part" => "index.md",
        "Table of Content" => "toc.md",
        "Tree Functionality" => "test1.md",
        "Likelihood" => "test2.md",
        "Parser" => "parser.md"]
   )

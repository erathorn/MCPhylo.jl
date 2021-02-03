include("../MCPhylo.jl")
using .MCPhylo

trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

Plots.plot(trees[1], treetype=:fan, msc=:black, mc=:green)

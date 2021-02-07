include("../MCPhylo.jl")
using .MCPhylo

trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

plot1 = Plots.plot(trees[1], showtips=true, msc=:green, mc=:yellow, lc=:white, bg=:black, tipfont=(7, :white))

include("../MCPhylo.jl")
using .MCPhylo

trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

plot(trees[1], :fan, msc=:green)
plot("./doc/Tree/Drav_mytrees_1.nwk", msc=:purple)

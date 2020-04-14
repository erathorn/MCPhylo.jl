
include("../MCPhylo.jl")
include("../Tree/Tree_Traversal.jl")

# COMPILER DOESN'T UNDERSTAND US UNLESS WE WILL DO THIS ^^^^

# OKAY, SO the minimal thing like this will work!
F = parsing_the_newick("(A,B,E)F;",nothing,0)
println(F.children)

# BUT THIS ALREADY RETURNS BULLSHIT

F = parsing_the_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F; ",nothing,0)
println(F.children)

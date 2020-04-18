
#include("../MCPhylo.jl")
#include("ParseNewick.jl")
#include("../Tree/Node_Type.jl")
include("../Tree/Tree_Traversal.jl")
#include("../Tree/Tree_Basics.jl")


function runthis(newick::String,current_node::Any)
while length(newick)>1
    new_tmp, F = parsing_the_newick(string(newick),current_node,0)
    println(F.children)
    newick = new_tmp
    current_node = F
end
end
runthis("(A,(C)E)F; ",nothing)
# BUT THIS ALREADY RETURNS BULLSHIT
#
# F = parsing_the_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F; ",nothing,0)
# println(F.children)

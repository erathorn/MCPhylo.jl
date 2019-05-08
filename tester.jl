#=
tester:
- Julia version: 1.0.1
- Author: Johannes
- Date: 2019-05-07
=#


include("./my_tree.jl")
using .my_tree

this_tree = my_tree.create_tree_from_leaves([1,2,3])

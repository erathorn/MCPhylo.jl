#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#


include("./Tree_module.jl")
using .Tree_module
using DataFrames

this_tree = Tree_module.create_tree_from_leaves([1,2,3,4,5])

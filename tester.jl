#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#


include("my_module.jl")
using .M


this_tree = M.create_tree_from_leaves([1,2,3,4,5])

M.to_df(this_tree)

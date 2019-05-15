#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#


include("Tree_Module.jl")
using .Tree_Module
using DataFrames

this_tree = Tree_Module.create_tree_from_leaves([1,2,3,4,5])

my_df = Tree_Module.Converter.to_df(this_tree)

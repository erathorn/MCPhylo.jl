#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#


include("./Tree/Tree_Module.jl")
using .Tree_Module
using DataFrames

this_tree = Tree_Module.create_tree_from_leaves([1 ,2 ,3, 4,5, 6, 7, 8])
Tree_Module.NNI!(this_tree)
my_df = Tree_Module.Converter.to_df(this_tree)
st = Tree_Module.Converter.from_df(my_df)


for sym in [:foo, :bar, :baz]
        @eval function $(Symbol(string("test_$sym")))(n::Int64)
                for i in 1:n
                    println($sym)
                end
        end
end

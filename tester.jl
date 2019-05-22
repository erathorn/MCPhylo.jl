#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#


include("./Tree/Tree_Module.jl")
include("./Substitution/SubstitutionMat.jl")
include("./Likelhood/LikelihoodCalculator.jl")
include("./Sampler/SamplerFunctions.jl")
using .Tree_Module
using .SubstitutionMat
using .LikelihoodCalculator
using .SamplerFunctions
using DataFrames

this_tree = Tree_Module.create_tree_from_leaves([1 ,2 ,3, 4,5, 6, 7, 8])
Tree_Module.NNI!(this_tree)
my_df = Tree_Module.Converter.to_df(this_tree)
st = Tree_Module.Converter.from_df(my_df)
a = [0.3,0.2,0.5]
SamplerFunctions.SampleVectorAdjusted!(a,1, 0.1, 0.0, 1.0, 1.0)

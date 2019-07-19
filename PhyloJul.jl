module PhyloJul

using Markdown
using Random
using Distributions
using DataFrames
using LinearAlgebra
using SpecialFunctions


include("./Tree/Tree_Basics.jl")
include("./Tree/Converter.jl")
include("./Tree/Tree_moves.jl")

include("./Parser/ParseNexus.jl")

include("./Substitution/SubstitutionMat.jl")

include("./Utils/SIMD_Mat.jl")

include("./Likelihood/LikelihoodCalculator.jl")
include("./Likelihood/Prior.jl")





end # module

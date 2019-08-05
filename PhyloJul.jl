module PhyloJul

using Markdown
using Mamba
using Random
using Distributions
using DataFrames
using LinearAlgebra
using SpecialFunctions



include("./Utils/Typing.jl")
include("./Tree/Tree_Basics.jl")
include("./Tree/Converter.jl")
include("./Tree/Tree_moves.jl")
include("./Tree/Tree_Matrix.jl")

include("./Parser/ParseNexus.jl")

include("./Sampler/ProbPathHMC.jl")
include("./Sampler/PhyloHMC_Functions.jl")


include("./Substitution/SubstitutionMat.jl")

include("./Utils/SIMD_Mat.jl")

include("./Likelihood/LikelihoodCalculator.jl")
include("./Likelihood/Prior.jl")





end # module

module PhyloJul

using Markdown
using Mamba
using Random
using Distributions
using DataFrames
using LinearAlgebra
using SpecialFunctions


abstract type ProbPathTune <: SamplerTune end
abstract type TreeVariate end


mutable struct TreeStochastic <: ScalarVariate
    value::TreeVariate
    symbol::Symbol
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::Distribution
end

mutable struct ProbPathHMCTune <: ProbPathTune
    n_leap::Int64
    stepsz::Float64

    ProbPathHMCTune() = new()
end # mutable struct

mutable struct Node <: TreeVariate
    name::String
    data::Array{Float64,2}
    child::Vector{Node}
    nchild::Int64
    root::Bool
    inc_length::Float64
    binary::String
    num::Int64

    Node() = new("",zeros(Float64,(1,2)), Node[], 0, true, 0.0, "0", 0)

    function Node(n::String, d::Array{Float64,2}, c::Vector{Node}, n_c::Int64, r::Bool, inc::Float64, b::String, num::Int64)
        mn = Node()
        mn.name = n
        mn.data = d
        mn.child = c
        mn.nchild = n_c
        mn.root = r
        mn.inc_length = inc
        mn.binary = b
        mn.num = num
        return mn
    end
end # struct Node



struct SamplerVariateP
    value::Node
    tune::ProbPathHMCTune

    function SamplerVariateP(v, t::ProbPathHMCTune)
        new(v, t)
    end
end


mutable struct SamplerP{T<:SamplerTune}
    params::Symbol
    eval::Function
    tune::T
    #targets::Vector{Symbol}
end


include("./Utils/Typing.jl")
include("./Tree/Tree_Basics.jl")
include("./Tree/Converter.jl")
include("./Tree/Tree_moves.jl")
include("./Tree/Tree_Matrix.jl")

include("./Parser/ParseNexus.jl")

include("./Sampler/ProbPathHMC.jl")

include("./Substitution/SubstitutionMat.jl")

include("./Utils/SIMD_Mat.jl")

include("./Likelihood/LikelihoodCalculator.jl")
include("./Likelihood/Prior.jl")

export size



end # module

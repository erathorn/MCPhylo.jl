module MCPhylo

using Reexport
@reexport using Distributions

#################### Imports ####################

using DelimitedFiles
using SpecialFunctions
using Serialization
using Distributed
using Printf: @sprintf
using LinearAlgebra
using Calculus: gradient
using Showoff: showoff
using Markdown
using DataFrames
using Random
using CSV

import Base: Matrix, names, summary
import Compose: Context, context, cm, gridstack, inch, MeasureOrNumber, mm, pt, px
import LinearAlgebra: cholesky, dot
import Statistics: cor
import Distributions:
       ## Generic Types
       Continuous, ContinuousUnivariateDistribution, Distribution,
       MatrixDistribution, Multivariate, MultivariateDistribution, PDiagMat,
       PDMat, ScalMat, Truncated, Univariate, UnivariateDistribution,
       ValueSupport,
       ## ContinuousUnivariateDistribution Types
       Arcsine, Beta, BetaPrime, Biweight, Cauchy, Chi, Chisq, Cosine,
       Epanechnikov, Erlang, Exponential, FDist, Frechet, Gamma, Gumbel,
       InverseGamma, InverseGaussian, Kolmogorov, KSDist, KSOneSided, Laplace,
       Levy, Logistic, LogNormal, NoncentralBeta, NoncentralChisq,
       NoncentralF, NoncentralT, Normal, NormalCanon, Pareto, Rayleigh,
       SymTriangularDist, TDist, TriangularDist, Triweight, Uniform, VonMises,
       Weibull,
       ## DiscreteUnivariateDistribution Types
       Bernoulli, Binomial, Categorical, DiscreteUniform, Geometric,
       Hypergeometric, NegativeBinomial, NoncentralHypergeometric, Pareto,
       PoissonBinomial, Skellam,
       ## MultivariateDistribution Types
       Dirichlet, Multinomial, MvNormal, MvNormalCanon, MvTDist,
       VonMisesFisher,
       ## MatrixDistribution Types
       InverseWishart, Wishart,
       ## Methods
       cdf, dim, gradlogpdf, insupport, isprobvec, logpdf, logpdf!, maximum,
       minimum, pdf, quantile, rand, sample!, support, length
import Gadfly: draw, Geom, Guide, Layer, layer, PDF, PGF, Plot, plot, PNG, PS,
       render, Scale, SVG, Theme
using LightGraphs: DiGraph, add_edge!, outneighbors,
       topological_sort_by_dfs, vertices
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats

include("distributions/pdmats2.jl")
using .PDMats2



"""
    Node

This data type holds the basic Node structure. The type T is used to specify the type of the data
stored in the node.

* If `nchild` is `0` the Node is a leaf node.
* If `root` is `False` the Node is a child of another node.
* `inc_length` specifies the length of the incomming branch.
* `binary` specifies the path from the root to the Node. `1` and `0` represent left and right turns respectively.
"""

mutable struct Node
    name::String
    data::Array{Float64,2}
    mother::Union{Node, Missing}
    lchild::Union{Node, Missing}
    rchild::Union{Node, Missing}
    nchild::Int64
    root::Bool
    inc_length::Float64
    binary::String
    num::Int64

    #Node() = new("",zeros(Float64,(1,2)), nothing, nothing, 0, true, 0.0, "0", 0)
    Node() = new("noname")
    function Node(n::String, d::Array{Float64,2}, m::Union{Node, Missing},c1::Union{Node, Missing},c2::Union{Node, Missing} , n_c::Int64, r::Bool, inc::Float64, b::String, num::Int64)
        mn = Node()
        mn.name = n
        mn.data = d
        mn.mother = m
        mn.lchild = c1
        mn.rchild = c2
        mn.nchild = n_c
        mn.root = r
        mn.inc_length = inc
        mn.binary = b
        mn.num = num
        return mn
    end
end # struct Node

#################### Types ####################

ElementOrVector{T} = Union{T, Vector{T}}


#################### Variate Types ####################

abstract type ScalarVariate <: Real end
abstract type ArrayVariate{N} <: DenseArray{Float64, N} end
abstract type TreeVariate <: Any end

const AbstractVariate = Union{ScalarVariate, ArrayVariate, TreeVariate}
const VectorVariate = ArrayVariate{1}
const MatrixVariate = ArrayVariate{2}



#################### Distribution Types ####################

const DistributionStruct = Union{Distribution,
                                 Array{UnivariateDistribution},
                                 Array{MultivariateDistribution}}


#################### Dependent Types ####################

mutable struct ScalarLogical <: ScalarVariate
  value::Float64
  symbol::Symbol
  monitor::Vector{Int}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
end

mutable struct ArrayLogical{N} <: ArrayVariate{N}
  value::Array{Float64, N}
  symbol::Symbol
  monitor::Vector{Int}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
end

mutable struct TreeLogical <: TreeVariate
  value::Node
  symbol::Symbol
  monitor::Vector{Int}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
end

mutable struct ScalarStochastic <: ScalarVariate
  value::Float64
  symbol::Symbol
  monitor::Vector{Int}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
  distr::UnivariateDistribution
end

mutable struct ArrayStochastic{N} <: ArrayVariate{N}
  value::Array{Float64, N}
  symbol::Symbol
  monitor::Vector{Int}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
  distr::DistributionStruct
end

mutable struct TreeStochastic <: TreeVariate
    value::Node
    symbol::Symbol
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::DistributionStruct
end

const AbstractLogical = Union{ScalarLogical, ArrayLogical, TreeLogical}
const AbstractStochastic = Union{ScalarStochastic, ArrayStochastic, TreeStochastic}
const AbstractDependent = Union{AbstractLogical, AbstractStochastic}
#const AnyDependent = Union{AbstractDependent, TreeStochastic}

#################### Sampler Types ####################

mutable struct Sampler{T}
  params::Vector{Symbol}
  eval::Function
  tune::T
  targets::Vector{Symbol}
end


abstract type SamplerTune end

struct SamplerVariate{T<:SamplerTune} <: VectorVariate
  value::Union{Vector{Float64}, Vector{Node}}
  tune::T

  function SamplerVariate{T}(x::AbstractVector, tune::T) where T<:SamplerTune
    v = new{T}(x, tune)
    validate(v)
  end

  function SamplerVariate{T}(x::AbstractVector, pargs...; kargs...) where T<:SamplerTune
    if !isa(x[1], Node)
      value = convert(Vector{Float64}, x)
    else
      value = convert(Vector{Node}, x)
    end
    new{T}(value, T(value, pargs...; kargs...))
  end
end


#################### Model Types ####################

struct ModelGraph
  graph::DiGraph
  keys::Vector{Symbol}
end

struct ModelState
  value::Vector{Any}
  tune::Vector{Any}
end


mutable struct Model
  nodes::Dict{Symbol, Any}
  samplers::Vector{Sampler}
  states::Vector{ModelState}
  iter::Int
  burnin::Int
  hasinputs::Bool
  hasinits::Bool
  likelihood::Float64
end


#################### Chains Type ####################

abstract type AbstractChains end

struct Chains <: AbstractChains
  value::Array{Float64, 3}
  range::StepRange{Int, Int}
  names::Vector{AbstractString}
  chains::Vector{Int}
  trees::Array{AbstractString, 3}
end

struct ModelChains <: AbstractChains
  value::Array{Float64, 3}
  range::StepRange{Int, Int}
  names::Vector{AbstractString}
  chains::Vector{Int}
  trees::Array{AbstractString, 3}
  model::Model
end


#################### Includes ####################

include("progress.jl")
include("utils.jl")
include("variate.jl")

include("distributions/constructors.jl")
include("distributions/distributionstruct.jl")
include("distributions/extensions.jl")
include("distributions/pdmatdistribution.jl")
include("distributions/transformdistribution.jl")
include("distributions/Phylodist.jl")
include("distributions/Branchlengthdist.jl")

include("model/dependent.jl")
include("model/graph.jl")
include("model/initialization.jl")
include("model/mcmc.jl")
include("model/model.jl")
include("model/simulation.jl")

include("output/chains.jl")
include("output/chainsummary.jl")
include("output/discretediag.jl")
include("output/fileio.jl")
include("output/gelmandiag.jl")
include("output/gewekediag.jl")
include("output/heideldiag.jl")
include("output/mcse.jl")
include("output/modelchains.jl")
include("output/modelstats.jl")
include("output/rafterydiag.jl")
include("output/stats.jl")
include("output/plot.jl")

include("samplers/sampler.jl")

include("samplers/abc.jl")
include("samplers/amm.jl")
include("samplers/amwg.jl")
include("samplers/bhmc.jl")
include("samplers/bia.jl")
include("samplers/bmc3.jl")
include("samplers/bmg.jl")
include("samplers/dgs.jl")
include("samplers/hmc.jl")
include("samplers/mala.jl")
include("samplers/miss.jl")
include("samplers/nuts.jl")
include("samplers/rwm.jl")
include("samplers/rwmc.jl")
include("samplers/slice.jl")
include("samplers/slicesimplex.jl")


include("Tree/Tree_Basics.jl")
include("Tree/Converter.jl")
include("Tree/Tree_moves.jl")
#include("Tree/Tree_Matrix.jl")

include("Parser/Parser.jl")
include("Parser/ParseCSV.jl")
include("Parser/ParseNexus.jl")

include("Sampler/ProbPathHMC.jl")
include("Sampler/ProbPathHMC_Node.jl")
include("Sampler/PhyloHMC_Functions.jl")
include("Sampler/PhyloHMC_Functions_Node.jl")
include("Sampler/BranchLengthSlice.jl")


include("Substitution/SubstitutionMat.jl")

include("Utils/SIMD_Mat.jl")
include("Utils/FileIO.jl")

#include("Likelihood/LikelihoodCalculator_Matrix.jl")
include("Likelihood/LikelihoodCalculator_Node.jl")
include("Likelihood/Prior.jl")
#################### Exports ####################

export
  AbstractChains,
  AbstractDependent,
  AbstractLogical,
  AbstractStochastic,
  AbstractVariate,
  ArrayLogical,
  ArrayStochastic,
  ArrayVariate,
  TreeLogical,
  TreeVariate,
  Node,
  Chains,
  Logical,
  MatrixVariate,
  Model,
  ModelChains,
  Sampler,
  SamplerTune,
  SamplerVariate,
  ScalarLogical,
  ScalarStochastic,
  ScalarVariate,
  Stochastic,
  VectorVariate

export
  BDiagNormal,
  Flat,
  SymUniform,
  CompoundDirichlet,
  PhyloDist,
  MultivariateUniformTrunc,
  CompoundDirichletWrap

export
  autocor,
  changerate,
  cor,
  describe,
  dic,
  discretediag,
  draw,
  gelmandiag,
  gettune,
  gewekediag,
  gradlogpdf,
  gradlogpdf!,
  graph,
  graph2dot,
  heideldiag,
  hpd,
  insupport,
  invlink,
  invlogit,
  link,
  logit,
  logpdf,
  logpdf!,
  mcmc,
  mcse,
  plot,
  predict,
  quantile,
  rafterydiag,
  rand,
  readcoda,
  relist,
  relist!,
  sample!,
  setinits!,
  setinputs!,
  setmonitor!,
  setsamplers!,
  summarystats,
  unlist,
  update!

@static if VERSION >= v"1.0.0"
  export showall
end

export
  ABC,
  AMM, AMMVariate,
  AMWG, AMWGVariate,
  BHMC, BHMCVariate,
  BMC3, BMC3Variate,
  BMG, BMGVariate,
  DiscreteVariate,
  DGS, DGSVariate,
  HMC, HMCVariate,
  BIA, BIAVariate,
  MALA, MALAVariate,
  MISS,
  NUTS, NUTSVariate,
  RWM, RWMVariate,
  Slice, SliceMultivariate, SliceUnivariate,
  SliceSimplex, SliceSimplexVariate,
  ProbPathHMC,
  BranchSlice,
  RWMC, RWMCVariate

export
  make_tree_with_data,
  to_file

export
  cm,
  inch,
  mm,
  pt,
  px


#################### Deprecated ####################

include("deprecated.jl")


end # module
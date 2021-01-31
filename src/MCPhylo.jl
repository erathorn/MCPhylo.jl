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
using Plots
using StatsPlots
using Zygote
using FiniteDiff
using Showoff: showoff
using Markdown
using DataFrames
using Random
using CSV
using ChainRules


using CUDA
if has_cuda()
  using GPUArrays
else
  @warn "The Julia CUDA library is installed, but no CUDA device detected.
         Computation is performed without CUDA functionality."
end


import Base: Matrix, names, summary, iterate
import Base.Threads.@spawn
import Compose: Context, context, cm, gridstack, inch, MeasureOrNumber, mm, pt, px
import LinearAlgebra: cholesky, dot, BlasFloat
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
using LightGraphs: DiGraph, add_edge!, outneighbors,
       topological_sort_by_dfs, vertices
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats

include("distributions/pdmats2.jl")
using .PDMats2
include("Tree/Node_Type.jl") # We need this to get the Node type in

#################### Types ####################

ElementOrVector{T} = Union{T, Vector{T}}


#################### Variate Types ####################

abstract type ScalarVariate <: Real end
abstract type ArrayVariate{N} <: DenseArray{Float64, N} end
abstract type TreeVariate <: AbstractNode end

const AbstractVariate = Union{ScalarVariate, ArrayVariate, TreeVariate}
const NumericalVariate = Union{ScalarVariate, ArrayVariate}
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

mutable struct TreeLogical{T} <: TreeVariate where T<:GeneralNode
  value::T
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
  value::Union{CuArray{Float64, N}, Array{Float64, N}}
  symbol::Symbol
  monitor::Vector{Int}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
  distr::DistributionStruct
end


mutable struct TreeStochastic{T} <: TreeVariate where T<: GeneralNode
    value::T
    symbol::Symbol
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::DistributionStruct
end

const AbstractLogical = Union{ScalarLogical, ArrayLogical}
const AbstractStochastic = Union{ScalarStochastic, ArrayStochastic}
const AbstractTreeStochastic = Union{TreeLogical, TreeStochastic}
const AbstractDependent = Union{AbstractLogical, AbstractStochastic, AbstractTreeStochastic}


#################### Sampler Types ####################

mutable struct Sampler{T}
  params::Vector{Symbol}
  eval::Function
  tune::T
  targets::Vector{Symbol}
end


abstract type SamplerTune end

struct SamplerVariate{T<:SamplerTune} <: VectorVariate
  value::Union{Vector{Float64}, Vector{S}} where S<:GeneralNode
  tune::T

  function SamplerVariate{T}(x::AbstractVector, tune::T) where T<:SamplerTune
    v = new{T}(x, tune)
    validate(v)
  end

  function SamplerVariate{T}(x::AbstractVector, pargs...; kargs...) where T<:SamplerTune
    if !isa(x[1], GeneralNode)
      value = convert(Vector{Float64}, x)
    else
      mt = typeof(x[1])
      value = convert(Vector{mt}, x)
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
  value::Vector
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
  names::Vector{S} where S <: AbstractString
  chains::Vector{Int}
  trees::Array{S, 3} where S <: AbstractString
  moves::Array{Int, 1}
end

struct ModelChains <: AbstractChains
  value::Array{Float64, 3}
  range::StepRange{Int, Int}
  names::Vector{S} where S <: AbstractString
  chains::Vector{Int}
  model::Model
  trees::Array{S, 3} where S <: AbstractString
  moves::Array{Int, 1}
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
include("distributions/TreeConstraints.jl")

include("model/dependent.jl")
include("model/dependent_tree.jl")
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
include("Tree/Tree_Distance.jl")
include("Tree/Tree_Traversal.jl")
include("Tree/Tree_Search.jl")
include("Tree/Tree_Legacy.jl")
include("Tree/Tree_Clustering.jl")
include("Tree/Tree_Ladderizing.jl")
include("Tree/Tree_Pruning.jl")
include("Tree/Tree_Consensus.jl")


include("Parser/Parser.jl")
include("Parser/ParseCSV.jl")
include("Parser/ParseNexus.jl")
include("Parser/ParseNewick.jl")
include("Sampler/PNUTS.jl")

include("Utils/FileIO.jl")

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
  TreeStochastic,
  GeneralNode,
  AbstractNode,
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
  exponentialBL

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
  PNUTS, PNUTSVariate

export
  make_tree_with_data,
  make_tree_with_data_cu,
  to_file, drop_samples,
  to_df,
  to_covariance,
  to_distance_matrix,
  node_distance,
  get_path,
  cut_samples,
  NNI!, NNI,
  SPR!, SPR,
  slide!, slide,
  swing!, swing,
  RF, randomize!,
  BHV_bounds,
  get_branchlength_vector,
  set_branchlength_vector!,
  post_order,
  pre_order,
  level_order,
  add_child!,
  remove_child!,
  delete_node!,
  insert_node!,
  find_lca,
  find_by_binary,
  find_by_name,
  find_root,
  create_tree_from_leaves,
  create_tree_from_leaves_cu,
  newick,
  tree_length,
  tree_height,
  path_length,
  get_sister,
  get_leaves,
  neighbor_joining,
  upgma,
  prune_tree!, prune_tree,
  ladderize_tree!, ladderize_tree,
  majority_consensus_tree


export
  cm,
  inch,
  mm,
  pt,
  px

#################### Deprecated ####################
println("test")

# bemp = parsing_newick_string("(((alak1253,cuaa1241,((((gola1254,kont1244,bahn1262)),(reng1252,seda1262)),(hala1252,jehh1245)),((east2333,bulo1242),sree1244),(nyah1249,lave1249,lave1248)),(paco1243,(ngeq1245,(lowe1395,ongg1239)),east1236,((east2332,kata1264),sooo1254)),((warj1242,(khas1269,pnar1238)),((manm1238,(west2396,para1301)),((ruch1235,shwe1236),rian1261))),(cent1989,nort2684),(((yuan1241,khmu1256)),mall1246),(mang1378,buga1247),(monn1252,nyah1250),((khar1287,(pare1266,sora1254)),sant1410),(carn1240,cent1990),((male1282,rucc1239),viet1252)),idum1241);")
#
# bendo = find_by_name(bemp,"idum1241")
# demp = reroot(bemp,"idum1241")


# tree,df = fill_in_the_blanks("C:/Users/Adham/Desktop/testdata.paps.nex", "(((alak1253,cuaa1241,((((gola1254,kont1244,bahn1262)),(reng1252,seda1262)),(hala1252,jehh1245)),((east2333,bulo1242),sree1244),(nyah1249,lave1249,lave1248)),(paco1243,(ngeq1245,(lowe1395,ongg1239)),east1236,((east2332,kata1264),sooo1254)),((warj1242,(khas1269,pnar1238)),((manm1238,(west2396,para1301)),((ruch1235,shwe1236),rian1261))),(cent1989,nort2684),(((yuan1241,khmu1256)),mall1246),(mang1378,buga1247),(monn1252,nyah1250),((khar1287,(pare1266,sora1254)),sant1410),(carn1240,cent1990),((male1282,rucc1239),viet1252)),idum1241);")
# node = find_by_name(tree,"khmu1256")
# println(newick(tree))
#
# println(df)


for file in cd(readdir, "C:/Users/Adham/Desktop/thesis_data")
  if occursin("con.tre", file)


# open("C:/Users/Adham/Desktop/shorterdata.paps.nex", "r") do file
#       global content = readlines(file)
#     end
#
#
# while true
#   line = popfirst!(content)
#   if line == "matrix"
#     break
#   end
# end
# tempdf = create_nexusdf(content)
# println(tempdf[1][1], "\t", tempdf[2][1])

# bla = make_tree_with_data("C:/Users/Adham/Desktop/shorterdata.paps.nex";replace_missing=false)
# da_dataframe = bla[2]
# da_tree = parsing_newick_string("(((alak1253,cuaa1241,((((gola1254,kont1244,bahn1262)),(reng1252,seda1262)),(hala1252,jehh1245)),((east2333,bulo1242),sree1244),(nyah1249,lave1249,lave1248)),(paco1243,(ngeq1245,(lowe1395,ongg1239)),east1236,((east2332,kata1264),sooo1254)),((warj1242,(khas1269,pnar1238)),((manm1238,(west2396,para1301)),((ruch1235,shwe1236),rian1261))),(cent1989,nort2684),(((yuan1241,khmu1256)),mall1246),(mang1378,buga1247),(monn1252,nyah1250),((khar1287,(pare1266,sora1254)),sant1410),(carn1240,cent1990),((male1282,rucc1239),viet1252)),idum1241);")
#
#
# for x in eachindex(da_dataframe[2,1,:])
#   # println("THIS IS NODE NUMBER ",x)
#   # println("\n\n")
#   for val in eachindex(da_dataframe[2,:,x])
#     if da_dataframe[2,val,x] == 3.0
#       # println("UNKNOWN VALUE DETECTED")
#       # println("index of unknown value is: ",[2,val,x])
#       global index_to_change = [2,val,x]
#     end #if
#   end #for
#   # println("NEW NODE \n \n STARTS NOW")
# end #for
# n_char = length(da_dataframe[2,:,1])
# dim1,dim2,dim3 = index_to_change
# global df_copy = []
# for node in eachindex(da_dataframe[2,1,:])
#   if da_dataframe[2,dim2,node] == da_dataframe[dim1,dim2,dim3] && [2,dim2,node] != [dim1,dim2,dim3]
#     # print("NODE NUMBER ")
#     # print(node)
#     # print(" HAS ? IN SAME PLACE AS NODE NUMBER ")
#     # println(dim3)
#   else
#     # push!(df_copy,da_dataframe[:,:,node])
#   end #if
# end #for
# # println(da_dataframe)
include("deprecated.jl")

end # module

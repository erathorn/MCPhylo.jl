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
@reexport using Plots
using Plots.PlotMeasures
@reexport using StatsPlots
@reexport using MCPhyloTree
import MCPhyloTree: internal_external
using RecipesBase
import RecipesBase.plot
using StatsBase
import StatsBase: predict
using Zygote
using FiniteDiff
using ForwardDiff
using ReverseDiff
using Showoff: showoff
using Markdown
using Random
using DataStructures
using ProgressMeter
using Bijectors
using LogExpFunctions
using LoopVectorization
using Base.Threads
using PDMats
using IfElse
using DiffResults
import Base: Matrix, names, summary


import Statistics: cor
import Distributions:
    #        ## Generic Types
    Continuous,
    ContinuousUnivariateDistribution,
    Distribution,
    MatrixDistribution,
    Multivariate,
    MultivariateDistribution,
    PDiagMat,
    PDMat,
    ScalMat,
    Truncated,
    Univariate,
    UnivariateDistribution,
    ValueSupport,
    # Functions
    cdf,
    dim,
    gradlogpdf,
    insupport,
    isprobvec,
    logpdf,
    logpdf!,
    maximum,
    minimum,
    pdf,
    quantile,
    rand,
    sample!,
    support,
    length

using LightGraphs: DiGraph, add_edge!, outneighbors, topological_sort_by_dfs, vertices
import StatsBase:
    autocor,
    autocov,
    countmap,
    counts,
    describe,
    predict,
    quantile,
    sample,
    sem,
    summarystats

import Bijectors: link

#################### Types ####################

ElementOrVector{T} = Union{T,Vector{T}}


#################### Variate Types ####################

abstract type AbstractVariate end


#################### Distribution Types ####################

const DistributionStruct =
    Union{Distribution,Array{UnivariateDistribution},Array{MultivariateDistribution}}


#################### Dependent Types ####################

mutable struct Logical{D<:Union{Real,AbstractArray{T,N} where {T<:Real,N},GeneralNode}, F<:Function} <:
    AbstractVariate
    value::D
    symbol::Symbol
    monitor::Vector{Int}
    eval::F
    sources::Vector{Symbol}
    targets::Vector{Symbol}
end



mutable struct Stochastic{
    D<:Union{Real,AbstractArray{T,N} where {T<:Real,N},GeneralNode},F<:Function, DI<:DistributionStruct
} <: AbstractVariate
    value::D
    symbol::Symbol
    monitor::Vector{Int}
    eval::F
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::DI
    lpdf::Float64
end

const ScalarVariate = Union{Stochastic{<:Real, <:Function, <:DistributionStruct},Logical{<:Real, <:Function}}
const VectorVariate =
    Union{Stochastic{<:AbstractArray{<:Real,1}, <:Function, <:DistributionStruct},
    Logical{<:AbstractArray{<:Real,1}, <:Function}}
const MatrixVariate =
    Union{Stochastic{<:AbstractArray{<:Real,2}, <:Function, <:DistributionStruct},
    Logical{<:AbstractArray{<:Real,2}, <:Function}}
const TreeVariate = Union{Logical{<:GeneralNode, <:Function},Stochastic{<:GeneralNode,<:Function,<:DistributionStruct}}

# General Union
const ArrayVariate = Union{
    Stochastic{<:AbstractArray{<:Real,N} where N, <:Function, <:DistributionStruct},
    Logical{<:AbstractArray{<:Real,N} where {N}, <:Function},
}


const AbstractLogical = Union{Logical{<:Real, <:Function},Logical{<:AbstractArray{<:Real,N} where {N}, <:Function}}
const AbstractStochastic =
    Union{Stochastic{<:Real, <:Function, <:DistributionStruct},Stochastic{<:AbstractArray{<:Real,N} where {N}, <:Function, <:DistributionStruct}}
const AbstractDependent = Union{AbstractLogical,AbstractStochastic,TreeVariate}

const ConstraintDict{S} = Dict{
    Symbol,
    Union{Vector{Vector{S}},Vector{Tuple{Vector{S},Vector{S}}}} where S<:AbstractString,
}

abstract type classic_nuts end
abstract type riemann_nuts end
abstract type fwd end
abstract type zyg end
abstract type rev end
abstract type fin end
abstract type provided end
abstract type nograd end

const NUTS_Form = Union{classic_nuts, riemann_nuts}
const GradType = Union{fwd, zyg, fin, rev, nograd, provided}

#################### Sampler Types ####################

abstract type SamplerTune end
abstract type GradSampler{G} <:SamplerTune end
mutable struct Sampler{
    T<:SamplerTune,
    R<:AbstractArray{S} where {S<:Union{Real,GeneralNode}},
} <: AbstractVariate
    value::R
    params::Vector{Symbol}
    tune::T
    targets::Vector{Symbol}
    transform::Bool
end




Base.length(s::Sampler) = length(s.value)

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
    nodes::Dict{Symbol,Any}
    samplers::Vector{Sampler}
    states::Vector{ModelState}
    iter::Int
    burnin::Int
    hasinputs::Bool
    hasinits::Bool
    likelihood::Float64
end


############## Additional Structs ################
struct SimulationParameters
    burnin::Int64
    thin::Int64
    chains::Int64
    verbose::Bool
    trees::Bool
    asdsf::Bool
    freq::Int64
    min_splits::Float64
end

mutable struct ConvergenceStorage
    splitsQueue::Vector{Accumulator{BitVector,Int64}}
    splitsQueues::Vector{Vector{Accumulator{BitVector,Int64}}}
    n_tree_gens::Int64
    ASDSF_vals::Vector{Vector{Float64}}
end

function ConvergenceStorage(splitsQueue, splitsQueues)
    ConvergenceStorage(splitsQueue, splitsQueues, 0, [Float64[] for i in splitsQueues])
end

#################### Chains Type ####################

abstract type AbstractChains end

struct Chains <: AbstractChains
    value::Array{Float64,3}
    range::StepRange{Int,Int}
    names::Vector{AbstractString}
    chains::Vector{Int}
    trees::Array{AbstractString,3}
    tree_names::Vector{AbstractString}
end

struct ModelChains <: AbstractChains
    value::Array{Float64,3}
    range::StepRange{Int,Int}
    names::Vector{AbstractString}
    chains::Vector{Int}
    model::Model
    trees::Array{AbstractString,3}
    tree_names::Vector{AbstractString}
    stats::Array{Float64,2}
    stat_names::Vector{AbstractString}
    sim_params::SimulationParameters
    conv_storage::Union{Nothing,ConvergenceStorage}
    samplers::Vector{Vector{Sampler}}
end

#################### Includes ####################

include("customerrors.jl")
include("utils.jl")
include("variate.jl")

include("distributions/TreeDistribution.jl")
include("distributions/distributionstruct.jl")
include("distributions/extensions.jl")
include("distributions/pdmatdistribution.jl")
include("Likelihood/SubstitutionModels.jl")
include("distributions/Phylodist.jl")
include("distributions/TreeConstraints.jl")

include("model/dependent.jl")
include("model/dependent_tree.jl")
include("model/graph.jl")
include("model/initialization.jl")
include("model/mcmc.jl")
include("model/model.jl")
include("model/simulation.jl")
include("model/simulation_statistics.jl")

include("output/asdsf.jl")
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

include("samplers/tree_hamiltonian/states.jl")
include("samplers/tree_hamiltonian/refraction.jl")
include("samplers/tree_hamiltonian/auxilliary.jl")

include("samplers/abc.jl")

include("samplers/dmh.jl")
include("samplers/empirical.jl")
include("samplers/hmc.jl")
include("samplers/miss.jl")
include("samplers/nuts.jl")
include("samplers/nuts_classical.jl")
include("samplers/nuts_riemannian.jl")
include("samplers/pnuts.jl")
include("samplers/pphmc.jl")
include("samplers/rwm.jl")
include("samplers/slice.jl")
include("samplers/slicesimplex.jl")

include("Parser/Parser.jl")
include("Parser/ParseCSV.jl")
include("Parser/ParseNexus.jl")

include("Likelihood/Prior.jl")
include("Likelihood/Rates.jl")
include("Likelihood/VectorizedFunctions.jl")
include("Likelihood/LikelihoodCalculator_Node.jl")
#################### Exports ####################

export AbstractChains,
    AbstractDependent,
    AbstractLogical,
    AbstractStochastic,
    AbstractVariate,
    ArrayVariate,
    ConstraintDict,
    ConvergenceStorage,
    TreeVariate,
    GeneralNode,
    AbstractNode,
    Node,
    Chains,
    FileSyntaxError,
    Logical,
    MatrixVariate,
    Model,
    ModelChains,
    Sampler,
    SamplerTune,
    ScalarVariate,
    SimulationParameters,
    Stochastic,
    TreeDistribution,
    VectorVariate

export BDiagNormal,
    Flat, SymUniform, CompoundDirichlet, PhyloDist, MultiplePhyloDist, exponentialBL, UniformBranchLength, UniformTopology, UniformConstrained, SymBinary

export ASDSF,
    PNUTS_monitor,
    autocor,
    changerate,
    contourplot,
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
    link,
    logpdf,
    logpdf!,
    mcmc,
    mcse,
    plot,
    plot_asdsf,
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

export ABC,
    HMC,
    HMCVariate,
    MISS,
    NUTS,
    NUTSVariate,
    RWM,
    RWMVariate,
    Slice,
    SliceMultivariate,
    SliceUnivariate,
    SliceSimplex,
    SliceSimplexVariate,
    DMH,
    DMHVariate,
    PNUTS,
    PNUTSVariate,
    PPHMC,
    PPHMCVariate,
    Empirical,
    EmpiricalVariate

export make_tree_with_data,
    to_file,
    drop_samples,
    generate_constraints,
    generate_constraints!,
    topological,
    discrete_gamma_rates,
    Restriction,
    JC,
    GTR,
    freeK,
    classic_nuts, riemann_nuts,
    fwd, zyg, rev, provided, fin, nograd


export cm, inch, mm, pt, px

end # module

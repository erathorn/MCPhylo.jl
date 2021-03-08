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
include("Utils/confusionutility.jl")

include("Likelihood/LikelihoodCalculator_Node.jl")
include("Likelihood/Prior.jl")
include("Likelihood/Rates.jl")
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

"""
language family data
"""
aa_datafile = "C:/Users/Adham/Desktop/thesis_data/aa/data-aa-58-200.paps.nex"
aa_newicktree = "((monn1252:0.005954094,nyah1250:0.005954094):0.03040002,(buga1247:0.02182118,mang1378:0.02182118):0.01453293,((((khmu1256:0.002660002,yuan1241:0.002660002):0.01733502,(mall1246:0.01505151,idum1241:0.01505151):0.004943511):0.01119221,(((west2396:0.004767014,para1301:0.004767014):0.01085023,(((ruch1235:0.005374165,shwe1236:0.005374165):0.004329879,rian1261:0.009704043):0.002997801,manm1238:0.01270184):0.002915401):0.01027385,(warj1242:0.00758499,(pnar1238:0.002663731,khas1269:0.002663731):0.004921259):0.0183061):0.005296147):0.002101967,((khar1287:0.01545796,(pare1266:0.00668051,sora1254:0.00668051):0.008777447):0.003207969,sant1410:0.01866593):0.01462328):0.003064912,(((((((irrr1240:0.003190416,lowe1395:0.003190416):0.001193639,ongg1239:0.004384055):0.006064092,(ngeq1245:0.008668157,(east1236:0.006642721,paco1243:0.006642721):0.002025435):0.001779991):0.002256727,((east2332:0.002390148,sooo1254:0.002390148):0.002312446,kata1264:0.004702594):0.00800228):0.01218627,((viet1252:0.005846834,rucc1239:0.005846834):0.004803253,male1282:0.01065009):0.01424106):0.005365286,((((lave1249:0.005756634,(nyah1249:0.00321849,lave1248:0.00321849):0.002538144):0.00838779,alak1253:0.01414442):0.003814294,(((((bahn1262:8.903838E-4,gola1254:8.903838E-4):8.072645E-4,kont1244:0.001697648):0.008566601,((hala1252:0.003999205,jehh1245:0.003999205):0.003625497,(seda1262:0.004283164,reng1252:0.004283164):0.003341539):0.002639547):0.001918164,cuaa1241:0.01218241):0.002647151,(sree1244:0.007865651,(east2333:0.005595922,bulo1242:0.005595922):0.002269729):0.006963913):0.003129154):0.008230523,(cent1989:0.001505817,nort2684:0.001505817):0.02468342):0.00406719):0.002062014,(carn1240:0.01200406,cent1990:0.01200406):0.02031439):0.004035672);"
aa_rates = "C:/Users/Adham/Desktop/thesis_data/aa/data-aa-58-200.paps.nex.pstat"

ie_datafile = "C:/Users/Adham/Desktop/thesis_data/ie/data-ie-42-208.paps.nex"
ie_newicktree = "((anci1242:0.01590758,mode1248:0.01590758):0.03272224,(((((russ1263:0.007192853,(bela1254:0.003942959,ukra1253:0.003942959):0.003249894):0.004145972,poli1260:0.01133883):0.004001035,((slov1269:0.004131983,czec1258:0.004131983):0.006620253,(uppe1395:0.003020021,lowe1385:0.003020021):0.007732216):0.004587624):0.005389967,(((bulg1262:0.004143415,mace1250:0.004143415):0.005634082,serb1264:0.009777497):0.005920406,slov1268:0.0156979):0.005031924):0.002615481,chur1257:0.02334531):0.02528452,(osse1243:0.03555512,((urdu1245:0.01035955,mait1250:0.01035955):0.008004383,mara1378:0.01836394,oriy1255:0.01836394):0.01719118):0.0130747,(((((stan1290:0.008506885,ital1282:0.008506885):0.004745095,((port1283:0.004505445,stan1288:0.004505445):0.004662041,stan1289:0.009167486):0.004084494):0.01000715,lati1261:0.02325913):0.01554813,oldi1245:0.03880725):0.004257193,(stan1293:0.02985418,(((((faro1244:0.005404019,(icel1247:0.002727279,oldn1244:0.002727279):0.002676739):0.004629725,norw1258:0.01003374):0.002952381,((swed1254:0.007534174,elfd1234:0.007534174):0.002192903,laum1238:0.009727077):0.003259048):0.002615087,dani1285:0.01560121):0.007833912,(stan1295:0.01049717,dutc1256:0.01049717):0.01293795):0.00641906):0.01321026):0.005565377);"
ie_rates = "C:/Users/Adham/Desktop/thesis_data/ie/data-ie-42-208.paps.nex.pstat"
ie_rates2 = "C:/Users/Adham/Desktop/thesis_data/ie/data-ie-42-208.paps.nex (1).pstat"
ie_newicktree2 = "(stan1290:0.006838919,ital1282:0.003780526,(((port1283:0.00463633,stan1288:0.002768069):0.003549467,stan1289:0.004301185):0.002008977,((((stan1293:0.008391929,(((((faro1244:0.002620913,(icel1247:0.0023488,oldn1244:0.002225082):0.002241154):0.004585207,norw1258:0.001382665):0.001467979,((swed1254:0.002669744,elfd1234:0.01286755):0.001158366,laum1238:0.002654463):0.001373923):0.001506777,dani1285:0.001639329):0.009633214,(stan1295:0.00571627,dutc1256:0.005710192):0.005887555):0.007536698):0.03402907,(((((((russ1263:0.005320513,(bela1254:0.002751365,ukra1253:0.002842417):0.001663942):0.00234324,poli1260:0.00613651):0.00187831,((slov1269:0.002001475,czec1258:0.002354197):0.002751759,(uppe1395:7.096634E-4,lowe1385:0.002864164):0.005675873):0.00194154):0.003247953,(((bulg1262:0.003302828,mace1250:0.002545882):0.003566056,serb1264:0.002956689):0.002778568,slov1268:0.00381578):0.001267526):0.00144318,chur1257:0.002287584):0.04235327,(osse1243:0.07755994,(((urdu1245:0.01474154,mait1250:0.006753825):0.006130628,mara1378:0.01554617):0.004295805,oriy1255:0.008663987):0.05490788):0.03230486):0.004275028,(anci1242:8.942463E-4,mode1248:0.02489026):0.04999762):0.006480824):0.004091097,oldi1245:0.0653026):0.02413651,lati1261:0.00198977):0.01098717):0.001540279);"

st_datafile = "C:/Users/Adham/Desktop/thesis_data/st/data-st-64-110.paps.nex"
st_newicktree = "(nort2732:0.1083415,dhim1246:0.1083415,idum1241:0.1083415,(rawa1265:0.01404208,drun1238:0.01404208):0.09429945,(((chon1286:0.01640531,mong1332:0.01640531):0.04507972,((tang1336:0.04009198,(zeme1240:0.01465064,lian1251:0.01465064):0.02544134):0.01084329,loth1237:0.05093527):0.01054977):0.01662227,(hmar1241:0.01636882,lush1249:0.01636882):0.06173849):0.03023422,((((((ersu1241:0.01801074,lisu1250:0.01801074):0.02142027,nort2722:0.03943101):0.009487086,((shix1238:0.02873905,namu1246:0.02873905):0.01015407,guiq1238:0.03889313):0.01002497):0.009849024,(horp1239:0.04076874,quey1238:0.04076874):0.01799838):0.01545851,((lich1241:0.00844172,naxi1245:0.00844172):0.03866265,(((zhab1238:0.01150843,pela1242:0.01150843):0.01812008,nucl1310:0.02962851,youl1235:0.02962851):0.007799446,nusu1239:0.03742795):0.00967642):0.02712126):0.01783753,((((kham1282:0.009255027,(tibe1272:0.00640825,amdo1237:0.00640825):0.002846776):0.01509792,(sher1255:0.01370404,jire1238:0.01370404):0.01064891):0.04622216,dakp1242:0.07057511):0.01155904,(guru1261:0.06027297,kaik1246:0.06027297):0.02186117):0.009929014):0.01627837,((((mind1253:0.009466597,hakk1236:0.009466597):0.006040841,mand1415:0.01550744):0.009281254,oldc1244:0.02478869):0.05750113,(((sgaw1245:0.01487829,geba1237:0.01487829):0.008332487,(yint1235:0.01116863,(west2409:0.006903287,manu1255:0.006903287):0.004265344):0.01204214,(geko1235:0.01143928,yinb1236:0.01143928):0.0117715):0.00408721,paok1235:0.02729798):0.05499183):0.02605171,(sher1257:0.09242828,(((apat1240:0.014763,nyis1236:0.014763):0.006506275,bori1243:0.02126927):0.005749,pada1257:0.02701827):0.06541001,(kach1280:0.06180453,((dima1251:0.0130707,bodo1269:0.0130707):0.01131343,garo1247:0.02438413):0.0374204):0.03062375,(lepc1244:0.0735681,((kulu1253:0.02907011,khal1275:0.02907011):0.009650421,limb1266:0.03872053):0.03484757):0.01886018):0.01591325);"
st_rates = "C:/Users/Adham/Desktop/thesis_data/st/data-st-64-110.paps.nex.pstat"




"""
Used to test rerooting, parsing
"""
# udmurt = "((Nenets,Selkup)no,((Mansi,Khanty)no,((Mari,((Permyak,Zyrian)no,Udmert)no)no,(Hungarian,(Saami,((Estonian,Finnish)no,Mordvin)no)no)no)no)no)no;"
#
# udmert_tree = parsing_newick_string(udmurt)
# println(newick(udmert_tree))
# println("\n")
#
#
# rerooted_udmert = reroot(udmert_tree, "Udmert")
# println(newick(rerooted_udmert))





"""
code below generates confusion matrix + summary statistics from the three language families' files
"""
thingstotest = ((aa_datafile,aa_newicktree,aa_rates),(ie_datafile,ie_newicktree2,ie_rates2),(st_datafile,st_newicktree,st_rates))
final_results = []

for langdata in thingstotest
  datafile,newicktree,ratesfile = langdata
  langfam = split(datafile,"-")[2]
  println("CURRENT LANGUAGE: ", langfam)
  results,pairwise_results, tree = fill_in_the_blanks_test(datafile, ratesfile, newicktree)
  amnttrue = length(filter(a-> a == true, results))
  amntfalse = length(filter(a->a == false, results))




  confusion = generate_confusion_matrix(pairwise_results)
  accuracy,precision,fscore = get_summary_statistics(confusion)

  push!(final_results,(langfam,(accuracy,precision,fscore),confusion))
end #for

for x in final_results
  println("RESULTS FOR ", x[1])
  println("ACCURACY: ", x[2][1])
  println("PRECISION: ", x[2][2])
  println("FSCORE: ", x[2][3])
  println("RESULTING MATRIX: ")
  println(x[3][1,:])
  println(x[3][2,:])
end #for
println("[1,1] in matrix denotes true positives, [2,1] false positives, [2,2] true negatives, [1,2] false negatives")




include("deprecated.jl")

end # module

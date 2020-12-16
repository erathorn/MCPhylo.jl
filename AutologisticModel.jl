include("./src/MCPhylo.jl")
using .MCPhylo
import .MCPhylo: logcond
using Random
import Distributions: logpdf
Random.seed!(42)

using LightGraphs, Random
using LinearAlgebra, Distributions, StatsBase, StatsFuns

include("wals_features.jl")
include("neighbour_graphs.jl")
include("AutologistcDistr.jl")

lmat = create_linguistic_nmat(dc_langs)
nmat = create_nmat(dmat, 500)
nmat = Int64.(nmat)

ngraph = SimpleGraph(nmat)
lgraph = SimpleGraph(lmat)

ov_space = ov_concordant_sums(data_array, nmat)
ov_ling = ov_concordant_sums(data_array, lmat)
ov_uni = universality(data_array)

spacesums = cond_concordant_sums(data_array, ngraph)
lingsums = cond_concordant_sums(data_array, lgraph)

nvals = maximum(data_array[1,:])

my_data = Dict{Symbol, Any}(
  :df => data_array[1,:]
);

ov_spacesum = ov_space[1]
ov_lingsum = ov_ling[1]
ov_unisum = ov_uni[1,:]

cond_space = spacesums[1,:,:]
cond_ling = lingsums[1,:,:]

nlangs = length(data_array[1,:])
nfeat = 1

# model setup
#where is the data?
model =  Model(
    df = Stochastic(2, (linw, spaw, uniw) ->
        AutologisticDistr(cond_ling, ov_lingsum, linw, cond_space, ov_spacesum, spaw,
		ov_unisum, uniw, nvals, nlangs, nfeat), false, false),
	linw = Stochastic(1, ()-> Gamma(1.0, 1.0), true),
	spaw = Stochastic(1, ()-> Gamma(1.0, 1.0), true),
	uniw = Stochastic(2, ()-> Normal(0, 10), true)
     )

# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(
 :df => data_array,
 :lingsums => lingsums,
 :spacesums => spacesums,
 :ov_ling => ov_ling,
 :ov_space => ov_space,
 :ov_uni => ov_uni,
 :linw => rand(4),
 :spaw => rand(4),
 :uniw => rand(4,9))
 for i in 1:2
]

scheme = [MCPhylo.DMH(:linw, 3000, 1.0), #m and s
 		MCPhylo.DMH(:spaw, 3000, 1.0),
		MCPhylo.DMH(:uniw, 3000, 1.0)
           ]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 10000, #n generations
   burnin=100,thin=100, chains=2)

to_file(sim, "newerDMH")
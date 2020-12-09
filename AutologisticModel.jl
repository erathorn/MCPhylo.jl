include("./src/MCPhylo.jl")
using .MCPhylo
import .MCPhylo: logcond
using Random
import Distributions: logpdf
Random.seed!(42)

using LightGraphs, Random
using LinearAlgebra, Distributions, StatsBase

include("wals_features.jl")
include("neighbour_graphs.jl")
include("AutologistcDistr.jl")

lmat = create_linguistic_nmat(langs)
nmat = create_nmat(sample_dmat, 1000)
nmat = Int64.(nmat)

ngraph = SimpleGraph(nmat)
lgraph = SimpleGraph(lmat)

ov_space = ov_concordant_sums(data_array, nmat)
ov_ling = ov_concordant_sums(data_array, lmat)
ov_uni = universality(data_array)

spacesums = cond_concordant_sums(data_array, ngraph)
lingsums = cond_concordant_sums(data_array, lgraph)

nvals = extract_nvals(data_array)

my_data = Dict{Symbol, Any}(
  :df => data_array
);

nlangs = length(data_array[1,:])
nfeat = 4

# model setup
model =  Model(
    df = Stochastic(2, (linw, spaw, uniw) ->
        AutologisticDistr(lingsums, ov_ling, linw, spacesums, ov_space, spaw,
		ov_uni, uniw, nvals, nlangs, nfeat), false, false),
	linw = Stochastic(1, ()-> Gamma(1.0, 1.0), true),
	spaw = Stochastic(1, ()-> Gamma(1.0, 1.0), true),
	uniw = Stochastic(1, ()-> Normal(0, 10), true)
     )


# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(
 :df => data_array,
 :lingsums => lingsums,
 :spacesums => spacesums,
 :ov_ling => ov_ling,
 :ov_space => ov_space,
 :ov_uni => ov_uni,
 :linw => rand(183),
 :spaw => rand(183),
 :uniw => rand(183))
 for i in 1:10]

scheme = [MCPhylo.DMH(:linw, 250, 1.0), #m and s
 		MCPhylo.DMH(:spaw, 250, 1.0),
		MCPhylo.DMH(:uniw, 250, 1.0)
           ]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 10, #n generations
   burnin=2,thin=1, chains=10)

to_file(sim, "myDMH")

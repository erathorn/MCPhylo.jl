
include("./src/MCPhylo.jl")
using .MCPhylo
import .MCPhylo: logcond
using Random
import Distributions: logpdf
Random.seed!(42)

using LinearAlgebra, Distributions, StatsBase, StatsFuns

include("wals_features.jl")
include("neighbour_graphs.jl")
include("AutologistcDistr.jl")

lmat = Float64.(create_linguistic_nmat(d))
nmat = create_nmat(dmat, 1000)

spacesums = cond_concordant_sums(data_array, nmat)
lingsums = cond_concordant_sums(data_array, lmat)

nvals = transpose(maximum(data_array, dims=2))[1,:]
my_data = Dict{Symbol, Any}(
  :df => data_array
);

nlangs = length(data_array[1,:])
nfeat = 1
uni_inits = uni_probs(data_array)

model =  Model(
    df = Stochastic(2, (linw, spaw, uniw) ->
        AutologisticDistr(lingsums, linw, lmat, spacesums, spaw, nmat,
		uniw, nvals, nlangs, nfeat), false, false),
	linw = Stochastic(1, ()-> Gamma(1.0, 1.0), true),
	spaw = Stochastic(1, ()-> Gamma(1.0, 1.0), true),
	uniw = Stochastic(2, ()-> Normal(0, 10), true)
     )

# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(
 :df => data_array,
 :lingsums => lingsums,
 :spacesums => spacesums,
 :linw => rand(1),
 :spaw => rand(1),
 :uniw => uni_inits)
 for i in 1:2]

scheme = [MCPhylo.DMH(:linw, 2700, 1.0, false, link=log, inverselink=exp),
 		MCPhylo.DMH(:spaw, 2700, 1.0, false, link=log, inverselink=exp),
		MCPhylo.DMH(:uniw, 2700, 0.01, false)
           ]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 50000, burnin=5000,thin=1, chains=2)

to_file(sim, "DMH02")

include("./src/MCPhylo.jl")
using .MCPhylo
import .MCPhylo: logcond
using Random
import Distributions: logpdf
Random.seed!(42)

using LightGraphs, MetaGraphs, Random
using LinearAlgebra, Distributions, StatsBase

include("wals_features.jl")
include("neighbour_graphs.jl")

lmat = create_linguistic_nmat(langs)
nmat = create_nmat(sample_dmat, 1000)
nmat = Int64.(nmat)

ngraph = MetaGraph(nmat)
lgraph = MetaGraph(lmat)
set_k_vals!(datadict, ngraph)
set_k_vals!(datadict, lgraph)

# to retrieve the value for language 1, feature k:
get_prop(ngraph, 5, Symbol("81A"))

"""
function for retrieving all the neighbour sums for x's and k's and
returning them as an array of (nlanguages, nfeatures) such that
the sums array[i,k] returns the sum for the ith language and the
kth feature. (The names of the features should still be passed
as a vector of strings).
"""

function n_sum_all(X::Array{Union{Missing,Int64},2}, g::MetaGraph, features::Vector{String})
	nvals = maximum(skipmissing(X)) # problem: might get number of features if nfeatures > max value
	nfeatures, nlang = size(X)
	sums = Array{Float64,3}(undef, nfeatures, nlang, nvals)
	for feature_idx in 1:nfeatures
		sym = Symbol(features[feature_idx])
		ith_feature = X[feature_idx,:]
		# get language index and language feature value
		for (i, x) in enumerate(ith_feature)
			for k in 1:nvals
				sums[feature_idx, i, k] = neighbour_k_sum(X[feature_idx,:], g, i, k, sym)
			end
		end
	end
	return sums
end

spatial_sums = n_sum_all(data_array, ngraph, feature_list)
linguistic_sums = n_sum_all(data_array, lgraph, feature_list)

include("AutologistcDistr.jl")

data_array[ismissing.(data_array)] .= -100

my_data = Dict{Symbol, Any}(
  :df => data_array
);


# model setup
model =  Model(
    df = Stochastic(2, (linw, spaw, uniw) ->
        AutologisticDistr(linguistic_sums, linw, spatial_sums, spaw, uniw, 183, 4), false, false),
	linw = Stochastic(1, ()->Normal(), true),
	spaw = Stochastic(1, ()->Normal(), true),
	uniw = Stochastic(1, ()->Normal(), true)
     )


# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
 :df => data_array,
 :linguistic_sums => linguistic_sums,
 :spatial_sums => spatial_sums,
 :linw => rand(183),
 :spaw => rand(183),
 :uniw => rand(183)
 ),
]


scheme = [MCPhylo.DMH(:linw, 250, 1.0)
           ]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 5,
   burnin=1,thin=1, chains=1)

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
		max_val = maximum(ith_feature)
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
#sums[feature,language,value]

"""
size(wo_co_sum) => 1296, 7
size(nusy_co_sum) => 1296, 4
ov_con_sum -> 3d array with #feauters x #languages x #(max(feature_vals))+1
[1, 3,:end] => 7
[2, 3, :end] => 4

"""
# To do
function weighted_neighbour_sums(x::Array{Int64,1}, ling_graph::MetaGraph,
	spatial_graph::MetaGraph)
	return x
end

"""
cond_sum takes the neighbour sums calculated above; in the sampler,
we will loop through the array and normalise the probabilities.
Args:
- the data array X
- the language index l
- the feature index f
- the spatial sums
- the linguistic sums
- the simulated parameters v, h, u
"""

function cond_sum(X::Array{Union{Missing,Int64},2}, l::Int64, f::Int64, spatial_sums::Array{Float64,3},
	ling_sums::Array{Float64,3}, v_params, h_params, u_params)
	nvals = maximum(skipmissing(X[f,:]))
	prob_kth_val = Vector{Float64}()
	for k in 1:nvals
		p = v_params[l] * spatial_sums[f,l,k] + h_params[l] * ling_sums[f,l,k] + u_params[l]
		push!(prob_kth_val, p)
	end
	return exp.(prob_kth_val)
end

"""
Inner sampler for DMH - args:
	- X: the data vector,
	- v_params: array of generated parameters v',
	- h_params: array of generated parameters h',
	- u_params: array of generated parameters u',
	- spatial_sums, ling_sums: neighbour sums for language l, feature k,
	- m: the number of update steps.
"""

r = range(0.0001, 0.05, length=100) # some range
dummy_v = rand(r, length(1:183)) # dummy parameters
dummy_h = rand(r, length(1:183))
dummy_u = rand(r, length(1:183))

function inner_sampler(X::Array{Union{Missing,Int64},2}, v_params, h_params,
	u_params, spatial_sums::Array{Float64,3}, ling_sums::Array{Float64,3}, m::Int64)
	nfeatures, nlang = size(X)
	@assert m >= nlang
	samples = deepcopy(X)
	counter = 0
	while true
		random_language_idx = shuffle(1:nlang)
		random_feature_idx = shuffle(1:nfeatures)
		for i in random_language_idx
			for f in random_feature_idx
				@show f
				feature_vals = sort!(unique(skipmissing(X[f,:])))
				@show feature_vals
				if ismissing(X[i])
					missing_probs = cond_sum.(Ref(samples), i, f, Ref(spatial_sums),
					Ref(ling_sums), Ref(v_params), Ref(h_params), Ref(u_params))
					missing_probs ./= sum(missing_probs)
					new_x = sample(feature_vals, StatsBase.weights(missing_probs))
					samples[f, i] = new_x
				end
					probs = cond_sum.(Ref(samples), i, f, Ref(spatial_sums), Ref(ling_sums),
						Ref(v_params), Ref(h_params), Ref(u_params))
					probs ./= sum(probs)
					@show probs
					new_x = sample(feature_vals, StatsBase.weights(probs))
					samples[f, i] = new_x
				end
			end
			counter += 1
			if counter > m
				return samples
			end
		end
	end
end


# TEST
inner_sampler(data_array, dummy_v, dummy_h, dummy_u, spatial_sums, linguistic_sums, 183)

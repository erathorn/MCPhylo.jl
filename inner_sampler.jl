using LightGraphs, MetaGraphs, Random
using LinearAlgebra, Distributions, StatsBase


"""
Including the two files below loads sorted data and functions for
creating neighbourhood matrices, to enable you to test the code
on some data (in this case, the word order).
"""

include("linguistic_features.jl")
include("neighbour_graphs.jl")

wo_nmat = create_nmat(wo_dmat, 500)
wo_nmat = Int64.(wo_nmat)
wo_lmat = create_linguistic_nmat(wo_data)

"""
neighbour_sum takes args:
- x: vector of data (feature values),
- nmat: either linguistic or spatial adjacency matrix,
- k: feature value of interest, and
- l: language (data point) of interest.

NB: it throws a BoundsError when I try to iterate through x in the sampler,
probably an issue with the way I index x[i], but I don't have
a good solution (except using graphs, which are easier to work with
due to the neighbor() function).
"""

function neighbour_sum(x::Array{Int64,1}, nmat::Array{Int64,2}, k::Int64, l::Int64)
	concordant_pairs = 0.0
	neighbour_indices = findall(x -> x != 0, nmat[l,:])
    for i in neighbour_indices
        if x[i] == k
			concordant_pairs += 1
		end
    end
    return concordant_pairs
end

"""
The functions below make use of graphs instead of matrices.
This will enable us to more easily implement weighted graphs
in the future.
"""

function set_xvals!(g::MetaGraph, X::Array{Int64,1})
    for (i, x) in enumerate(X)
        set_prop!(g, i, :x, x)
    end
    return g
end

function create_neighbour_graph(nmat::Array{Int64,2}, X::Array{Int64,1})
	g = MetaGraph(nmat)
	set_xvals!(g, X)
	return g
end

wo_ngraph = create_neighbour_graph(wo_nmat, wo_values)
wo_lgraph = create_neighbour_graph(wo_lmat, wo_values)

"""
concordant neighbour sum: the sum of all neighbours of l that share feature k
(when l also has feature k).
- X: data vector,
- g: spatial or linguistic graph,
- l: the  language (vertex in the graph),
- k: the feature value (e.g. SOV, SVO...)

neighbour k sum: the sum of all neighbours of l that have feature k
(when v doesn't necessarily have feature k). Same args as above.

Uncertainty about whether we should calculate only the concordant sums,
or the overall k sums.
"""

function concordant_neighbour_sum(X::Array{Int64,1}, g::MetaGraph, l::Int64, k::Int64)
	k_pairs = 0.0
	concordant_pairs = []
	for n in neighbors(g, v)
        if get_prop(g, n, :x) == get_prop(g, v, :x)
			pair = (n, v)
			push!(concordant_pairs, pair)
		end
	end
	for (p1, p2) in concordant_pairs
		if get_prop(g, p1, :x) == k
			k_pairs += 1
		end
		if get_prop(g, p2, :x) == k
			k_pairs += 1
		end
	end
	return k_pairs
end

function neighbour_k_sum(X::Array{Int64,1}, g::MetaGraph, l::Int64, k::Int64)
	sum = 0.0
	for n in neighbors(g, l)
		if get_prop(g, l, :x) == k
			sum += 1
		end
	end
	return sum
end

r = range(0.0001, 0.05, length=100) # some range
dummy_v = rand(r, length(wo_values)) # dummy parameters
dummy_h = rand(r, length(wo_values))
dummy_u = rand(r, length(wo_values))

"""
function for retrieving all the neighbour sums for x's and k's and
returning them as an array of (nlanguages, nfeatures) such that
the sums array[i,k] returns the sum for the ith language and the
kth feature.
"""

function get_neighbour_sums(X::Array{Int64,1}, g::MetaGraph)
	feature_vals = unique(X)
	N = length(X)
	K = length(feature_vals)
	sums = Array{Float64,2}(undef, N, K)
	for (i, x) in enumerate(X)
		for k in feature_vals
			sums[i,k] = neighbour_k_sum(X, g, i, k)
		end
	end
	return sums
end

# To do
function weighted_neighbour_sums(x::Array{Int64,1}, ling_graph::MetaGraph,
	spatial_graph::MetaGraph)
	return x
end

spatial_sums = get_neighbour_sums(wo_values, wo_ngraph)
linguistic_sums = get_neighbour_sums(wo_values, wo_lgraph)

"""
cond_prob takes the neighbour sums calculated above; in the sampler,
we will loop through it and normalise the probabilities.
"""

function cond_prob(X::Array{Int64,1}, l::Int64, k::Int64, spatial_sums::Array{Float64,2},
	ling_sums::Array{Float64,2}, v_params, h_params, u_params)
	p = v_params[l] * spatial_sums[l,k] + h_params[l] * ling_sums[l,k] + u_params[l]
    return exp(p)
end

"""
Inner sampler for double MH - args:
	- X: the data vector,
	- v_params: array of generated parameters v',
	- h_params: array of generated parameters h',
	- u_params: array of generated parameters u',
	- spatial_sums, ling_sums: neighbour sums for language l, feature k,
	- m: the number of update steps.
"""

function inner_sampler(X::Array{Int64,1}, v_params, h_params,
	u_params, spatial_sums::Array{Float64,2}, ling_sums::Array{Float64,2}, m::Int64)
	feature_vals = sort!(unique(X))
	N = length(X)
	@assert m >= N
	samples = deepcopy(X)
	idx_list = Vector(1:N)
	counter = 0
	while true
		random_index_order = shuffle(1:N)
		for i in random_index_order
			if ismissing(X[i]) # no missing values yet, and this part needs improvement;
				# To do: add the neighbour and global majority methods
				throw("We should not end up here")
				missing_probs = cond_prob.(Ref(samples), i, feature_vals, Ref(spatial_sums),
					Ref(ling_sums), Ref(v_params), Ref(h_params), Ref(u_params))
				missing_probs ./= sum(missing_probs)
				new_x = sample(feature_vals, StatsBase.weights(missing_probs))
				samples[i] = new_x
			end
				probs = cond_prob.(Ref(samples), i, feature_vals, Ref(spatial_sums), Ref(ling_sums),
					Ref(v_params), Ref(h_params), Ref(u_params))
				probs ./= sum(probs)
				new_x = sample(feature_vals, StatsBase.weights(probs))
				samples[i] = new_x
			end
			counter += 1
			if counter > m
				return samples
			end
		end
	end
end


# TEST
inner_sampler(wo_values, dummy_v, dummy_h, dummy_u, spatial_sums, linguistic_sums, 1296)

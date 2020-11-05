using LightGraphs, MetaGraphs
using LinearAlgebra, Distributions, StatsBase
using Random
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
neighbour_k_sum: the sum of all neighbours of v that share feature k
(when v also has feature k).
- X: data vector,
- g: spatial or linguistic graph,
- v: the vertex (i.e. language),
- k: the feature value (e.g. SOV, SVO...)
"""

function neighbour_k_sum(X::Array{Int64,1}, g::MetaGraph, v::Int64, k::Int64)
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

"""
Conditional distribution: Probability that x[l][i] = k.
u[l][k] - the endogenous tendency of the feature to take value k.
NB: this is only the numerator of the conditional distribution;
in the samper, I will loop through and normalise the probabilities.
"""

function cond_prob(X::Array{Int64,1}, k::Int64, l::Int64, ling_graph::MetaGraph,
	spatial_graph::MetaGraph, v::Array{Float64,1}, h::Array{Float64,1}, u::Array{Float64,1})
	# l -> language
	# k -> feature_value
	V = neighbour_k_sum(X, ling_graph, l, k)
	H = neighbour_k_sum(X, spatial_graph, l, k)
    p = v[l] * V + h[l] * H + u[l]
    return exp(p)
end

"""
Inner sampler for double MH - args:
	- x: the data vector;
	- v: array of generated parameters v';
	- h: array of generated parameters h';
	- nmat, lmat: the adjacency matrices;
	- m: the number of update steps (or define based on some other var?)
"""

r = range(0.0001, 0.05, length=100) # some range
dummy_v = rand(r, length(wo_values)) # dummy parameters
dummy_h = rand(r, length(wo_values))
dummy_u = rand(r, length(wo_values))

"""
Randomise iteration:
- create an index list;
- pick an index randomly from 1:length(x); it must be one of the
values in x;
- subtract that index from the index list to avoid repeats.

Maybe there's a better way to randomise the order of iteration?
(Also, I'm pretty sure the samples list ends up being in some random order,
so it won't correspond to the order of the original list -
not sure how to fix this.)

Also, I'd like a more efficient way to normalise than creating two lists
of probabilities, if possible.
"""

idxlist = Vector(1:length(wo_values))
r = rand(idxlist)


function sample_x(x, v, h, u, spatial_graph, ling_graph, m)
	nlang, n = size(spatial_graph)
	feature_vals = sort!(unique(x))
	#samples = Array{Int64,1}()
	N = length(x)
	@assert m >= N
	samples = deepcopy(x)
	idx_list = Vector(1:N)
	counter = 0
	while true
	#while m != 0 # is this the best way to limit the number of iterations?
		random_index_order = shuffle(1:N)
		for i in random_index_order
			#i = rand(idx_list) # primitive randomisation - problem: returns a vector
			if ismissing(x[i]) # no missing values yet, and this part needs improvement;
				# To do: add the neighbour and global majority methods
				throw("We should not end up here")
				missing_probs = Array{Float64,1}()
				missing_probs_norm = Array{Float64,1}()
				for k in feature_vals
					p_unnorm = cond_prob(x, k, i, spatial_graph, ling_graph, dummy_v, dummy_h, dummy_u)
					push!(missing_probs, p_unnorm)
				end
				for p in missing_probs
					p_norm = p/sum(missing_probs)
					push!(missing_probs_norm, p_norm)
				end
				new_x = sample(feature_vals, StatsBase.weights(missing_probs_norm))
				push!(samples, new_x) # with random iteration, how do I make sure samples[i] corresponds to the right x[i]?
			else
				#probs = Array{Float64,1}()
				#probs_norm = Array{Float64,1}() # inefficient?


				probs = cond_prob.(Ref(samples), feature_vals, i, Ref(spatial_graph), Ref(ling_graph), Ref(dummy_v), Ref(dummy_h), Ref(dummy_u))
				probs ./= sum(probs) # same as original probs_norm

		#		normalizer = 0.0
		#		for k in feature_vals
			#		p_unnorm = cond_prob(x, k, i, spatial_graph, ling_graph, dummy_v, dummy_h, dummy_u)
			#		normalizer += p_unnorm
			#		push!(probs, p_unnorm)
				#end
				#probs ./= normalizer # same as original probs_norm

				#for p in probs
				#	p_norm = p/sum(probs)
				#	push!(probs_norm, p_norm)
				#end
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
sample_x(wo_values, dummy_v, dummy_h, dummy_u, wo_ngraph, wo_lgraph, 1296)

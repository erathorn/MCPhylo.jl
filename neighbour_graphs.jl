using Distances, LinearAlgebra, LightGraphs, MetaGraphs

function create_neighbour_list(dm::Array{Float64,2}, threshold::Int64, languages::Vector{String})
    neighbours = Vector()
    for (i, l1) in enumerate(languages)
        for (j, l2) in enumerate(languages)
            if (dm[i,j] <= threshold) && (dm[i,j] != 0)
                push!(neighbours, [l1,l2])
            end
        end
    end
    return neighbours
end

function create_nmat(dm::Array{Float64,2}, threshold::Int64)
    n, m = size(dm)
    nmat = zeros(n, m)
    idx_list = Vector(1:n)
    for (i, l1) in enumerate(idx_list)
        for (j, l2) in enumerate(idx_list)
            if (dm[i,j] <= threshold) && (dm[i,j] != 0)
                nmat[i,j] = 1
            end
        end
    end
    return nmat
end

function replace_diagonal!(nmat::Matrix)
    for i in diagind(nmat, 0)
        nmat[i] == 1 ? nmat[i] = 0 : nothing
    end
    return nmat
end

function create_linguistic_nmat(d::DataFrame)
    nmat = [Int(g1==g2) for g1 in d.Genus, g2 in d.Genus]
    replace_diagonal!(nmat)
    return nmat
end

# to do
function create_weighted_graph(x)
    return x
end

"""
The functions below find neighbouring pairs for vertex v
that are concordant for trait x.
"""

function concordant_sum(X::Array{Int64,1}, g::MetaGraph, v::Int64)
    concordant_sum = 0
    for n in neighbors(g, v)
#        if has_prop(mg, :x, v) - build this in to throw an error if prop not defined for g, v
        if get_prop(g, n, :x) == get_prop(g, v, :x)
            concordant_pairs += 1 # or another way to measure the pairs?
        end
    end
    return concordant_sum
end

"""
nmat_sum takes args:
- x: vector of data (feature values),
- nmat: either linguistic or spatial adjacency matrix,
- k: feature value of interest, and
- l: language (data point) of interest.

NB: it throws a BoundsError when I try to iterate through x in the sampler,
probably an issue with the way I index x[i], but I don't have
a good solution (except using graphs, which are easier to work with
due to the neighbor() function).
"""

function nmat_sum(x::Array{Int64,1}, nmat::Array{Int64,2}, k::Int64, l::Int64)
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
set_k_vals! sets the feature values for the MetaGraph. Args:
- the dictionary of data (symbols to feature value arrays);
- the MetaGraph.
"""

function set_k_vals!(data::Dict, g::MetaGraph)
	nfeatures = length(data)
	for key in keys(data)
		feature_vals = data[key]
    	for (v, x) in enumerate(feature_vals)
        	set_prop!(g, v, key, x)
		end
    end
    return g
end


function create_neighbour_graph(nmat::Array{Int64,2}, X::Array{Any,2})
	g = MetaGraph(nmat)
	set_k_vals!(g, X)
	return g
end

"""
concordant neighbour sum: the sum of all neighbours of l that share feature k
(when l also has feature k).
- X: data vector,
- g: spatial or linguistic graph,
- l: the  language (vertex in the graph),
- k: the feature value (e.g. SOV, SVO...)
- sym: the feature value name as a symbol, e.g. Symbol("81A")

neighbour k sum: the sum of all neighbours of l that have feature k
(when v doesn't necessarily have feature k). Same args as above.
"""

function concordant_neighbour_sum(X::Array{Int64,1}, g::MetaGraph, l::Int64, k::Int64, sym::Symbol)
	k_pairs = 0.0
	concordant_pairs = []
	for n in neighbors(g, v)
        if get_prop(g, n, sym) == get_prop(g, v, sym)
			pair = (n, v)
			push!(concordant_pairs, pair)
		end
	end
	for (p1, p2) in concordant_pairs
		if get_prop(g, p1, sym) == k
			k_pairs += 1
		end
		if get_prop(g, p2, sym) == k
			k_pairs += 1
		end
	end
	return k_pairs
end

function neighbour_k_sum(X::Array{Union{Int64,Missing},1}, g::MetaGraph, l::Int64, k::Int64, sym::Symbol)
	sum = 0.0
	for n in neighbors(g, l)
		val = get_prop(g, n, sym)
		if ismissing(val)
			continue
		elseif val == k
			sum += 1
		end
	end
	return sum
end

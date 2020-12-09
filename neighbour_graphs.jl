using Distances, DataFrames, LinearAlgebra, LightGraphs

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
count_concordant_edges calculates the total number of concordant edges in
an adjacency matrix. Args:
- x: vector of data (values for one feature),
- nmat: either linguistic or spatial adjacency matrix
"""

function count_concordant_edges(X::Array{Int64,1}, nmat::Array{Int64,2})
	sum = 0
	n = length(X)
	for i in 1:n, j in 1:n
	    if i < j && nmat[i,j] == 1
			if X[i] == X[j] && X[i] ≠ -10
				sum += 1
			end
		end
	end
	return sum
end

"""
function ov_concordant_sums: retrieves the overall number of concordant pairs
for each feature in a graph (regardless of value). Returns an array
of length(nfeatures) with one sum per feature.
"""

function ov_concordant_sums(X::Array{Int64,2}, nmat::Array{Int64,2})
	sums = Vector{Float64}()
	nfeatures, nlang = size(X)
	for feat in 1:nfeatures
		sum_feat = count_concordant_edges(X[feat,:], nmat)
		push!(sums, sum_feat)
	end
	return sums
end

"""
neighbour k sum: the sum of all neighbours of l that have feature value k
(when l doesn't necessarily have value k). Needed for the sums for
the conditional distribution.
"""

function neighbour_k_sum(X::Array{Int64,1}, g::SimpleGraph, l::Int64, k::Int64)
	sum = 0.0
	for n in neighbors(g, l)
		val = X[n]
		if val == -10
			sum += 0
		elseif val == k
			sum += 1
		end
	end
	return sum
end

"""
function cond_concordant_sums: retrieves all the neighbour sums for features
and values and returns them as an array of (nfeatures, nlangs, nvals)
such that the new array[i,l,k] returns the sum for the ith language and the
kth feature value.
"""

function cond_concordant_sums(X::Array{Int64,2}, g::SimpleGraph)
	nvals = maximum(X)
	nfeatures, nlang = size(X)
	sums = Array{Float64,3}(undef, nfeatures, nlang, nvals)
	for feature_idx in 1:nfeatures
		lang_vals = X[feature_idx,:]
		for (lang_idx, x) in enumerate(lang_vals)
			for k in 1:nvals
				sums[feature_idx, lang_idx, k] = neighbour_k_sum(X[feature_idx,:], g, lang_idx, k,)
			end
		end
	end
	return sums
end

"""
universality: gets the universality of each feature value per feature.
Returns an array of dimensions nfeatures, nvals (max val).
arr[feature, value] = the sum of all langs for that feature with that value.
"""

function universality(X::Array{Int64,2})
	nfeatures, nlangs = size(X)
	nvals = maximum(X)
	sums = zeros(nfeatures, nvals)
	for i in 1:nfeatures
		feat_values = X[i,:]
		nvals_x = maximum(feat_values)
		for k in 1:nvals_x
			sum_k = count(x -> x == k, feat_values)
			sums[i, k] = sum_k
		end
	end
	return sums
end

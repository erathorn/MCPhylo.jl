using Distances, DataFrames, LinearAlgebra, LightGraphs, SimpleWeightedGraphs

include("wals_features.jl")

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

function standardize_rows!(nmat::Array{Float64,2})
	for i in 1:size(nmat,1) #rows
		for j in 1:size(nmat,2) #cols
			if nmat[i,j] == 0
				continue
			else
				row = nmat[i,:,:]
				col = nmat[:,j,:]
				nmat[i,:,:] = row ./ (sum(col) * nmat[i,j])
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

function create_weights_matrix(dm::Array{Float64,2})
	n, m = size(dm)
	wmat = zeros(n, m)
	N = Vector(1:n)
	for (i, v1) in enumerate(N)
		for (j, v2) in enumerate(N)
			if dm[i,j] == 0
				wmat[i, j] = 0
			else
				wmat[i,j] = 1/dm[i,j]
			end
		end
	end
    return wmat
end

"""
count_concordant_edges calculates the total number of concordant edges in
an adjacency matrix. Args:
- x: vector of data (values for one feature),
- nmat: either linguistic or spatial adjacency matrix
"""

function count_concordant_edges(X::Array{N,1}, nmat::Array{Int64,2})::Float64 where N <: Real
	conc_sum = zero(Float64)
	n = length(X)
	@inbounds for i in 1:n, j in 1:n
	    if i < j && nmat[i,j] == 1
			if X[i] == X[j] && X[i] ≠ -10
				conc_sum += 1
			end
		end
	end
	return conc_sum
end

"""
function ov_concordant_sums: retrieves the overall number of concordant pairs
for each feature in a graph (regardless of value). Returns an array
of length(nfeatures) with one sum per feature.

function weighted_ov_sum does the same, but uses the sum of the edge weights.
"""

function ov_concordant_sums(X::Array{N,2}, nmat::Array{Int64,2})::Vector{Float64} where N <: Real

	nfeatures, nlang = size(X)
	sums = Vector{Float64}(undef, nfeatures)
	@inbounds for feat in 1:nfeatures
		#sum_feat = count_concordant_edges(X[feat,:], nmat)
		#push!(sums, sum_feat)
		sums[feat] = count_concordant_edges(X[feat,:], nmat)
	end
	return sums
end

function find_concordant_pairs(X::Array{Int64,1})
	nlang = length(X)
	pairs = Vector{Vector{Int64}}()
	for i in 2:size(X,1)
    	for j in 1:(i-1)
        	if X[i] == X[j] && X[i] ≠ -10
				push!(pairs, [i,j])
			end
		end
	end
	return pairs
end

function ov_weighted_sum(X::Array{Int64,2}, wmat::Array{Float64,2})
	nfeat, nlang = size(X)
	sums = zeros(nfeat)
	for f in 1:nfeat
		xvals = X[f,:]
		pairs = find_concordant_pairs(xvals)
		for (i,pair) in enumerate(pairs)
			sums[f] += wmat[pair[1], pair[2]]
		end
	end
	return sums
end

"""
find_neighbours retrieves the indices of all the neighbours for each point (lang).

"""

function find_neighbours(X::Array{Int64,1}, nmat::Array{Float64,2}, lang::Int64; weighted=false)
	neighbours = Vector{Int64}()
	for (i,n) in enumerate(nmat[lang,:,:])
		if n ≠ 0
			push!(neighbours, i)
		end
	end
	return neighbours
end

function neighbour_k_sum(X::Array{Int64,1}, nmat::Array{Float64,2}, lang::Int64, k::Int64)
	sum = 0.0
	neighbours = find_neighbours(X,nmat,lang)
	for j in neighbours
		w = nmat[lang,j]
		if X[lang] == k
			sum += w
		end
	end
	return sum
end

"""
function cond_concordant_sums: retrieves all the neighbour sums for features
and values and returns them as an array of (nfeatures, nlangs, nvals)
such that the new array[f,l,k] returns the sum for the ith language and the
kth feature value. This can be used for either the weighted or regular
binary neighbourhood matrix.
"""

function cond_concordant_sums(X::Array{Int64,2}, nmat::Array{Float64,2})
	nvals = maximum(X)
	nfeatures, nlang = size(X)
	sums = zeros(nfeatures, nlang, nvals)
	for f in 1:nfeatures
		xvals = X[f,:]
		for (l, x) in enumerate(xvals)
			for k in 1:nvals
				sums[f, l, k] += neighbour_k_sum(xvals, nmat, l, k)
			end
		end
	end
	return sums
end

"""
universality: gets the universality sum of each feature value per feature.
Returns an array of dimensions nfeatures, nvals (max val).
arr[feature, value] = the sum of all langs for that feature with that value.
"""

function universality_sum(X::Array{Int64,2})
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

function uni_logprob(X::Array{Int64,2}) # for initial values of uni params
	nfeat, nlang = size(X)
	nvals = maximum(X)
	probs = zeros(nfeat, nvals)
	sums = universality_sum(X)
	for f in 1:nfeat
		feat_sums = sums[f,:]
		p = feat_sums./sum(feat_sums)
		probs[f,:] = p
	end
	return probs
end


mutable struct AutologisticDistr <: DiscreteMatrixDistribution
    ling_concordant::Array{Float64,3}
	ov_ling_concordant::Array{Float64,1}
    ling_params::Vector{Float64}
    spatial_concordant::Array{Float64,3}
	ov_spatial_concordant::Array{Float64,1}
    spatial_params::Vector{Float64} # length?
	ov_universality::Array{Float64,2}
    universality_params::Vector{Float64} # should params be vectors?
	nvals::Vector{Int64} # for each feature, the number of possible values
	nlangs::Int64
	nfeat::Int64
end # mutable struct

Base.size(d::AutologisticDistr) = (d.nfeat, d.nlangs)

"""
The unnormalised logpdf based on:
V[i,l,k] * v param[l] + H[i,l,k] * h param[l] + U[i,l,k] * u param[l]

It's necessary to loop through values as well as features since U is a 2D array
of dims (features, vals).
The parameters should be length(nlang).
"""

function logpdf(d::AutologisticDistr, X::Array{N,2}) where N <: Real
	res = 0
	maxval = maximum(d.nvals)
	for i in 1:d.nfeat # following logic in M & Y paper: V[i,l,k]
		for l in 1:d.nlangs
			res += d.ov_ling_concordant[i] * d.ling_params[l] +
			d.ov_spatial_concordant[i] * d.spatial_params[l]
			for k in 1:maxval
				# for some reason k comes out as a float unless converted?
				k = Int64(k)
				res += d.ov_universality[i,k] * d.universality_params[l]
			end
		end
	end
	return res
end

function logcond(d::AutologisticDistr, X::Array{N, 2}, l::Int64, f::Int64) where N <: Real
	# l -> language_index (is this being indexed in the sampler?)
	# k -> index of feature value
	nv = d.nvals[f] # -10, 1, 2, 4
	prob_kth_val = Vector{Float64}(undef, nv)
	for k in 1:nv
		p = d.spatial_params[l] * d.spatial_concordant[f,l,k] +
			d.ling_params[l] * d.ling_concordant[f,l,k] +
			d.universality_params[l]
		prob_kth_val[k] = exp(p)
	end
	return prob_kth_val
end

mutable struct ALRdistr <: DiscreteMatrixDistribution
    ling_concordant::Array{Float64,3}
	ov_ling_concordant::Array{Float64,1}
    ling_params::Vector{Float64}
    spatial_concordant::Array{Float64,3}
	ov_spatial_concordant::Array{Float64,1}
    spatial_params::Vector{Float64}
	ov_universality::Array{Float64,1}
    universality_params::Vector{Float64} # or: covars replace u
	covars::Vector{Float64}
	nlangs::Int64
	nfeat::Int64
end

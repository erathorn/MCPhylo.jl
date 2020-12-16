
mutable struct AutologisticDistr <: DiscreteMatrixDistribution
    ling_concordant::Array{Float64,3}
	ov_ling_concordant::Array{Float64,1}
    ling_params::Vector{Float64}
    spatial_concordant::Array{Float64,3}
	ov_spatial_concordant::Array{Float64,1}
    spatial_params::Vector{Float64} # length?
	ov_universality::Array{Float64,2}
    universality_params::Array{Float64,2} # should params be vectors?
	nvals::Vector{Int64} # for each feature, the number of possible values
	nlangs::Int64
	nfeat::Int64
end # mutable struct

Base.size(d::AutologisticDistr) = (d.nfeat, d.nlangs)

function logpdf(d::AutologisticDistr, X::Array{N,2}) where N <: Real
	res = 0
	maxval = maximum(d.nvals)
	for f in 1:d.nfeat
		res += d.ov_ling_concordant[f] * d.ling_params[f] +
		d.ov_spatial_concordant[f] * d.spatial_params[f]
		for k in 1:maxval
			res += d.ov_universality[f,k] * d.universality_params[f,k]
		end
	end
	return res
end

"""
P = exp(vertical_params[i] * vertical_sums[l,i,k] + ... + uni_params[i,k])
"""

using StatsFuns


function logcond(d::AutologisticDistr, X::Array{N, 2}, l::Int64, f::Int64) where N <: Real
	# l -> language_index (is this being indexed in the sampler?)
	# k -> index of feature value
	nv = d.nvals[f]
	prob_kth_val = Vector{Float64}(undef, nv)
	for k in 1:nv
		p = d.spatial_params[f] * d.spatial_concordant[f,l,k] +
			d.ling_params[f] * d.ling_concordant[f,l,k] +
			d.universality_params[f,k]
		prob_kth_val[k] = p
		end
	return softmax(prob_kth_val)
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

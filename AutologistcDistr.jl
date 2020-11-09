
mutable struct AutologisticDistr <: DiscreteMatrixDistribution
    ling_concordant::Array{Float64,3}
    ling_params::Vector{Float64}
    spatial_concordant::Array{Float64,3}
    spatial_params::Vector{Float64}
    universality_params::Vector{Float64}
	nlangs::Int64
	nfeat::Int64
end # mutable struct

Base.size(d::AutologisticDistr) = size(d.nlangs, d.nfeat)

"""
	This is the unnormalized logpdf
IT WILL NOT WORK
"""
function logpdf(d::AutologisticDistr, X::Array{N,2}) where N <: Real
	res = 0
	for i in 1:d.nfeat
		res += d.verticalWeights[i] * d.concordantVertical[i] + d.horizontalWeights[i] * d.concordantHorizontal[i] + d.universality[i]
	end
	res
end

function logcond(d::AutologisticDistr, X::Array{N, 2}, l::Int64, k::Int64) where N <: Real
	# l -> language_index
	# k -> index of feature value
	p = 0.0
	for i in 1:d.nfeat
		p += spatial_params[l] * spatial_concordant[l,i,k] + ling_params[l] * ling_concordant[l,i,k] + universality_params[l]
    end
    return exp(p)
end

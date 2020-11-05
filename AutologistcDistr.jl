
mutable struct AutologisticDistr <: DiscreteMatrixDistribution
    concordantVertical::Vector{Float64}
    verticalWeights::Vector{Float64}
    concordantHorizontal::Vector{Float64}
    horizontalWeights::Vector{Float64}
    universality::Vector{Float64}
	nlangs::Int64
	nfeat::Int64
end # mutable struct

Base.size(d::AutologisticDistr) = size(d.nlangs, d.nfeat)

"""
	This is the unnormalized logpdf
"""

function logpdf(d::AutologisticDistr, X::Array{N,2}) where N <: Real
	res = 0
	for i in 1:d.nfeat
		res += d.verticalWeights[i] * d.concordantVertical[i] + d.horizontalWeights[i] * d.concordantHorizontal[i] + d.universality[i]
	end
	res
end

function logcond(d::AutologisticsDistr, X::Array{N, 1}, l::Int64, k::Int64) where N <: Real
	# l -> language_index
	# k -> feature_value
	V = neighbour_k_sum(X, d.horizontalNeighbors, l, k)
	H = neighbour_k_sum(X, d.verticalNeigbors, l, k)
    p = d.verticalWeights[l] * V + d.horizontalWeights[l] * H + d.universality[l]
    return exp(p)
end

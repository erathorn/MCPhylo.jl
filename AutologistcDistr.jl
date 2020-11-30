
mutable struct AutologisticDistr <: DiscreteMatrixDistribution
    ling_concordant::Array{Float64,3}
    ling_params::Vector{Float64}
    spatial_concordant::Array{Float64,3}
    spatial_params::Vector{Float64}
    universality_params::Vector{Float64}
	nlangs::Int64
	nfeat::Int64
end # mutable struct

Base.size(d::AutologisticDistr) = (d.nfeat, d.nlangs)

"""
	This is the unnormalized logpdf
IT WILL NOT WORK
"""
function logpdf(d::AutologisticDistr, X::Array{N,2}) where N <: Real
	res = 0
	#for i in 1:d.nfeat
	#	res += d.verticalWeights[i] * d.concordantVertical[i] + d.horizontalWeights[i] * d.concordantHorizontal[i] + d.universality[i]
	#end
	res
end

function logcond(d::AutologisticDistr, X::Array{N, 2}, l::Int64, f::Int64) where N <: Real
	# l -> language_index
	# k -> index of feature value
	nvals = Int(maximum(skipmissing(X[f,:])))

	prob_kth_val = Vector{Float64}(undef, nvals)
	for k in 1:nvals
		p = d.spatial_params[l] * d.spatial_concordant[f,l,k] +
			d.ling_params[l] * d.ling_concordant[f,l,k] +
			d.universality_params[l]
		prob_kth_val[k] = exp(p)
	end
	return prob_kth_val
end

mutable struct ALRbinary <: DiscreteMatrixDistribution
    ling_sum::Array{Float64,2}
    ling_params::Vector{Float64}
    spatial_sum::Array{Float64,2} # sum of concordant minus discordant pairs
    spatial_params::Vector{Float64}
	covar_params::Vector{Float64}
	covars::Array{Float64,2} # n x p where p is the length of coefficient vector
	nlangs::Int64
	nfeat::Int64
end

function logpdf(d::ALRbinary, X::Array{N,2}) where N <: Real
	res = 0.0
	for f in 1:d.nfeat
		res += d.ling_sum[f] * d.ling_params[f] +
		d.spatial_sum[f] * d.spatial_params[f] +
		d.covars[] * d.covar_params[] # in the covars matrix, have something analogous to the sum
	end
	return exp(res)
end

function logcond(d::ALRbinary, lang::Int64)
	res = 0.0
	p = d.ling_params * d.ling_sum[i] +
		d.spatial_params * d.spatial_sum[i] + d.covar_params
	p = 2(p)
	return p
end

mutable struct AutologisticDistr <: DiscreteMatrixDistribution
	#cond_ling::Array{Float64,3}
#	ov_ling_concordant::Vector{Float64}
    ling_params::Vector{Float64}
	lingmat::Array{Float64,2}
    #cond_space::Array{Float64,3}
#	ov_spatial_concordant::Vector{Float64}
    spatial_params::Vector{Float64}
	spmat::Array{Float64,2}
#	ov_universality::Array{Float64,2}
    universality_params::Array{Float64,2}
	nvals::Vector{Int64}
	nlangs::Int64
	nfeat::Int64
end

Base.size(d::AutologisticDistr) = (d.nfeat, d.nlangs)

function logpdf(d::AutologisticDistr, X::Array{N,2}) where N <: Real
	res = zero(Float64)
	maxval = maximum(d.nvals)
	ling_cond = ov_concordant_sums(X, d.lingmat)
	sp_cond = ov_concordant_sums(X, d.spmat)
	#res = @. ling_cond[f] * d.ling_params[f] + sp_cond[f] * d.spatial_params[f]
	#println("ling: ", d.ling_params)
	@inbounds for f in 1:d.nfeat
		res += ling_cond[f] * d.ling_params[f] + sp_cond[f] * d.spatial_params[f]
		#println("---")
		Xf = X[f,:]
		@inbounds for k in 1:maxval
			uni_cond = count(x -> x == k, Xf)
			#println(uni_cond, " ", d.universality_params[f,k])
			res += uni_cond * d.universality_params[f,k]
		end
		#println("---")
	end
	#println("res $res")
	return res
end

"""
P = exp(vertical_params[i] * vertical_sums[l,i,k] + ... + uni_params[i,k])
"""

function logcond(d::AutologisticDistr, X::Array{N, 2}, l::Int64, f::Int64) where N <: Real
	nv = d.nvals[f]
	probs = d.universality_params[f, :]
	spsum = neighbour_sum(X, d.spmat, l, f, nv)
	lsum = neighbour_sum(X, d.lingmat, l, f, nv)
	@inbounds for k in 1:nv
		probs[k] += d.spatial_params[f] * spsum[k] +
					d.ling_params[f] * lsum[k]# +
	end
	return softmax(probs)
end

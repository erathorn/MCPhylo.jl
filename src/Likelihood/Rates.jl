"""
    function discrete_gamma_rates(α::T, β::S,k::Int64, method::Symbol=:mean; sig::Int64=1)::Array{Float64} where {T<:Real, S<:Real}

This function calculates the rates of the discretized gamma distribution.
It returns a vector of either the mean or the median weights in each category,
such that each category is of equal proportion.

See: Yang, 1994, Maximum likelihood phylogenetic estimation from DNA sequences
with variable rates over sites: Approximate methods. (https://doi.org/10.1007/BF00160154)
"""
function discrete_gamma_rates(α::T, β::S,k::Int64, method::Symbol=:mean; sig::Int64=1)::Array{Float64} where {T<:Real, S<:Real}
    meanvals = Array{Float64, 1}(undef, k)
    factor::Float64 = α/β*k

    if method==:median
        meanvals .= median_boundaries(α, β, k)
        t = sum(meanvals)
        meanvals .*= (factor / t)
    else
        meanvals .= mean_boundaries(α, β, k)
        @inbounds @simd for i in 1:k-1
            meanvals[i] = gamma_inc(α+1,meanvals[i]*β, sig)[1]
        end
        meanvals[k] = 1.0
        @inbounds for i in k:-1:2
            meanvals[i] -= meanvals[i-1]
            meanvals[i] *= factor
        end
        meanvals[1] *= factor
    end
    meanvals
end

#### helper functions for boundaries of rate categories ####

function mean_boundaries(α::T, β::S, k::Int64)::Array{Float64} where {T<:Real, S<:Real}
    ch = Chisq((2*α)/(2*β))
    boundaries = Array{Float64, 1}(undef, k)
    for i in 1:k-1
        boundaries[i] = quantile(ch, i/k)/(2*β)
    end
    boundaries
end

function median_boundaries(α::T, β::S, k::Int64)::Array{Float64} where {T<:Real, S<:Real}
    ch = Chisq((2*α)/(2*β))
    boundaries = Array{Float64, 1}(undef, k)
    for i in 1:k
        boundaries[i] = quantile(ch, ((i-1)*2+1)/(2*k))
    end
    return boundaries
end

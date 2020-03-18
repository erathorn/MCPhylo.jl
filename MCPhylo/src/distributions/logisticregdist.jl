
mutable struct BrownianPhylo <: DiscreteMatrixDistribution
        mu::Array{Float64,2}
        tree::Node
        Σ::Array{Float64}
        sigmai::Vector{Float64}
        P::Array{Float64,2}
        chars::Int64
        leaves::Int64
end

minimum(d::BrownianPhylo) = -Inf
maximum(d::BrownianPhylo) = +Inf
Base.size(d::BrownianPhylo) = (d.leaves, d.chars)

function logpdf(d::BrownianPhylo, x::AbstractArray{T, 2}) where T<:Real

    blv = get_branchlength_vector(d.tree)

    mlpdf(d.mu, d.tree, blv, d.sigmai, d.Σ, d.P, d.chars, x)

end

function mlpdf(mu::Array{Float64,2}, tree::Node, blv::Vector{T}, sigmai::Vector{Float64}, Σ::Array{Float64,2}, P::Array{Float64,2}, chars::Int64, data::AbstractArray{S, 2} where S) where T<:Real

    mycov::Array{T,2} = MCPhylo.to_covariance(tree, blv)

    rv::T = 0.0
    #fchars = Array{Function,2}(undef, size(data))
    #fchars .= x -> logcdf(Normal(),x)
    #fchars[data .== 1] .= x -> logccdf(Normal(), x)
    #si = sigmai .* Ref(Σ)

    #sis = si .* Ref(mycov)
    #noise = 1e-8.* Diagonal(ones(size(mu, 1)))
    ch = cholesky(mycov).L
    #si[i].*
    @inbounds @simd for i in 1:chars
        #ch = cholesky(sis[i]).L

        rv += sum(logpdf.(Bernoulli.(invlogit.(ch*P[:,i]+mu[:,i])), data[:,i]))


        #rv += sum(map((f,x)->f(x),fchars[:,i], nullstellen))::T
    end

    rv
end

function gradlogpdf(d::BrownianPhylo, x::AbstractArray{T, 2}) where T<:Real

    blv = get_branchlength_vector(d.tree)

    f(y) =mlpdf(d.mu, d.tree, y, d.sigmai, d.Σ, d.P, d.chars, x)

    r2 = DiffResults.GradientResult(blv)
    res = ForwardDiff.gradient!(r2, f, deepcopy(blv))


    DiffResults.value(res), DiffResults.gradient(res)
end

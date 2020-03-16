
mutable struct BrownianPhylo <: DiscreteMatrixDistribution
        mu::Array{Float64}
        tree::Node
        sigmai::Array{Float64}
        #P::Array{Float64}
        chars::Int64
        leaves::Int64
end

minimum(d::BrownianPhylo) = -Inf
maximum(d::BrownianPhylo) = +Inf
Base.size(d::BrownianPhylo) = (d.leaves, d.chars)

function logpdf(d::BrownianPhylo, x::AbstractArray{T, 2}) where T<:Real

    blv = get_branchlength_vector(d.tree)

    mlpdf(d.mu, d.tree, blv, d.sigmai, d.chars, x)

end

function mlpdf(mu::Vector{Float64}, tree::Node, blv::Vector{T}, sigmai::Vector{Float64}, chars::Int64, data::AbstractArray{S, 2} where S) where T<:Real

    nor = Normal()
    cov = MCPhylo.to_covariance(tree, blv)

    fchars = Array{Function,2}(undef, size(data))
    fchars .= x -> logcdf(nor,x)
    fchars[data .== 1] .= x -> logccdf(nor, x)
    cov += 1.0e-8 .* Diagonal(ones(size(cov, 1)))
    rv = 0.0
    mm = ones(size(data,1))
    #vs = [-2.33414, -0.741964, 0.741964, -2.33414]

    @simd for i in 1:chars
        ch = cholesky(sigmai[i]^2 .* cov).L
        nullstellen = ch \ -(mm .= mu[i])
        rv += sum(map((f,x)->f(x), fchars[:,i], nullstellen))
        #rv += sum(logpdf.(Bernoulli.(invlogit.(ch*P[:,i]+mu)), @view x[:,i]))
    end

    rv
end

function gradlogpdf(d::BrownianPhylo, x::AbstractArray{T, 2}) where T<:Real

    blv = get_branchlength_vector(d.tree)

    f(y) =mlpdf(d.mu, d.tree, y, d.sigmai, d.chars, x)

    r2 = DiffResults.GradientResult(blv)
    res = ForwardDiff.gradient!(r2, f, deepcopy(blv))


    DiffResults.value(res), DiffResults.gradient(res)
end

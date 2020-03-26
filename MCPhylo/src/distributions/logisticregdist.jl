
mutable struct BrownianPhylo <: DiscreteMatrixDistribution
        mu::Array{Float64,2}
        tree::Node{T,A,B,I} where {T<: Real, A<: AbstractArray, B<:AbstractArray, I<:Integer}
        sigmai::Vector{Float64}
        P::Array{Float64,2}
        scaler::Float64
        chars::Int64
        leaves::Int64
end

minimum(d::BrownianPhylo) = -Inf
maximum(d::BrownianPhylo) = +Inf
Base.size(d::BrownianPhylo) = (d.leaves, d.chars)

function logpdf(d::BrownianPhylo, x::AbstractArray{T, 2})::Float64 where T<:Real

    blv::Vector{Float64} = get_branchlength_vector(d.tree)

    mlpdf(d.mu, d.tree, blv, d.sigmai, d.P, d.scaler, d.chars, x)

end

function my_bernoulli_logpdf(θ::T, x::S)::T where {S <: Real, T<: Real}
    x == 1 ? θ : log(1-exp(θ))
end

function mlpdf(mu::Array{Float64,2}, tree::Node, blv::Vector{T}, sigmai::Vector{Float64},
               P::Array{Float64,2}, scaler, chars::Int64, data::Array{Float64,2})::T where T<:Real

    mycov::Array{T,2} = MCPhylo.to_covariance(tree, blv)
    rv::T = 0.0
    rv_a = zeros(T,chars)
    @inbounds @simd for i in 1:chars
        mat = (cholesky(sigmai[i].*mycov).L*P[:,i]+mu[:,i])
        rv_a[i] = sum(my_bernoulli_logpdf.(loginvlogit.(mat, Ref(scaler)), data[:,i]))
    end
    sum(rv_a)
end

function gradlogpdf(d::BrownianPhylo, x::AbstractArray{T, 2})::Tuple{Float64, Vector{Float64}} where T<:Real

    blv = get_branchlength_vector(d.tree)

    f(y) =mlpdf(d.mu, d.tree, y, d.sigmai, d.P,d.scaler, d.chars, x)

    r2 = DiffResults.GradientResult(blv)
    res = ForwardDiff.gradient!(r2, f, blv)

    grad = DiffResults.gradient(res)

    DiffResults.value(res), grad
end

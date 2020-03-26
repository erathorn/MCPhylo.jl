
mutable struct BrownianPhylo <: DiscreteMatrixDistribution
        mu::Array{Float64,2}
        tree::Node
        sigmai::Vector{Float64}
        P::Array{Float64,2}
        scaler::Float64
        chars::Int64
        leaves::Int64
end

minimum(d::BrownianPhylo) = -Inf
maximum(d::BrownianPhylo) = +Inf
Base.size(d::BrownianPhylo) = (d.leaves, d.chars)

function logpdf(d::BrownianPhylo, x::AbstractArray{T, 2}) where T<:Real

    blv = get_branchlength_vector(d.tree)

    mlpdf(d.mu, d.tree, blv, d.sigmai, d.P, d.scaler, d.chars, x)

end

function my_bernoulli_logpdf(θ::T, x::S)::T where {S <: Real, T<: Real}
    x == 1 ? θ : log(1-exp(θ))
end

function mlpdf(mu::Array{Float64,2}, tree::Node, blv::Vector{T}, sigmai::Vector{Float64},
      P::Array{Float64,2}, scaler, chars::Int64, data::AbstractArray{S, 2} where S) where T<:Real

    mycov::Array{T,2} = MCPhylo.to_covariance(tree, blv)
    rv::T = 0.0
    rv_a = zeros(T,chars)
    for i in 1:chars

        ch = cholesky(sigmai[i].*mycov).L
        mat = (ch*P[:,i]+mu[:,i])
        myvalue=sum(my_bernoulli_logpdf.(loginvlogit.(mat, Ref(scaler)), data[:,i]))
        #myvalue = sum(logpdf.(Bernoulli.(invlogit.(mat, Ref(scaler))), data[:,i]))

        rv_a[i] = myvalue

    end
    #println(rv_a)
    sum(rv_a)
end

function gradlogpdf(d::BrownianPhylo, x::AbstractArray{T, 2}) where T<:Real

    blv = get_branchlength_vector(d.tree)
    #println("blv ",blv)
    f(y) =mlpdf(d.mu, d.tree, y, d.sigmai, d.P,d.scaler, d.chars, x)
    #cfg1 = ForwardDiff.GradientConfig(f, blv, ForwardDiff.Chunk{1}());
    r2 = DiffResults.GradientResult(blv)
    res = ForwardDiff.gradient!(r2, f, blv)
    #println("blv ",blv)
    #println("grad ",DiffResults.gradient(res))
    grad = DiffResults.gradient(res)
    #grad_fd = Calculus.gradient(f, blv)
    #f(blv), grad_fd
    #println(isapprox(grad_fd , grad))
    #grad[isnan.(grad)] .= -Inf
    #if any(isnan.(grad))
        #println(blv)
    #    println(Calculus.gradient(f, blv))
        #grad .= -Inf
        #return -Inf, grad
        #throw("no way")

    #end
    DiffResults.value(res), grad#DiffResults.gradient(res)
end

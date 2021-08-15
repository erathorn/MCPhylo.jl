function internal_logpdf(d::CompoundDirichlet, b_lens::Array{Float64},
                         int_leave_map::Vector{Int64})

    blen_int = 0.0
    blen_leave = 0.0
    blen_int_log = 0.0
    blen_leave_log = 0.0
    nterm = 0.0


    @views for i in eachindex(int_leave_map)
        if int_leave_map[i] === 1
            @inbounds blen_int += b_lens[i]
            @inbounds blen_int_log += log(b_lens[i])
        else
            @inbounds blen_leave += b_lens[i]
            @inbounds blen_leave_log += log(b_lens[i])
            nterm += 1
        end
    end

    t_l = blen_int+blen_leave
    n_int = nterm-3.0

    first = (d.alpha*log(d.beta))-log(gamma(d.alpha)) - (t_l*d.beta)
    second = -log(gamma(d.a))-log(gamma(d.c))+log(gamma(d.a+d.c))
    third = blen_leave_log*(d.a-1.0) + blen_int_log*(d.a*d.c-1.0)
    fourth = (d.alpha-d.a*nterm-d.a*d.c*n_int)*log(t_l)

    r2 = first + second + third + fourth

    return r2

end

function gradlogpdf(d::CompoundDirichlet, x::FNode)
    int_ext = internal_external(x)
    blv = get_branchlength_vector(x)
    f(y) =  internal_logpdf(d, y, int_ext)
    r = Zygote.pullback(f, blv)
    return r[1],r[2](1.0)[1]
end

function gradlogpdf(d::exponentialBL, x::FNode)
    bl = get_branchlength_vector(x)
    g(y) = sum(logpdf.(Exponential(d.scale), y))
    r = Zygote.pullback(g, bl)
    r[1],r[2](1.0)[1]
end

function gradlogpdf(t::Union{UniformConstrained, UniformTopology, UniformBranchLength}, 
                    x::FNode)::Tuple{Float64, Vector{Float64}}

    blv = get_branchlength_vector(x)
    0.0, zeros(length(blv))
end

function logpdf(d::CompoundDirichlet, x::FNode)
    internal_logpdf(d, get_branchlength_vector(x), internal_external(x))
end

function logpdf(t::Union{UniformConstrained, UniformTopology, UniformBranchLength}, x::FNode)
    0.0
end

function logpdf(ex::exponentialBL, x::FNode)
    _logpdf(ex, x)
end

function _logpdf(d::exponentialBL, x::FNode)
    bl = get_branchlength_vector(x)
    sum(logpdf.(Exponential(d.scale), bl))
end

function logpdf_sub(d::CompoundDirichlet, x::FNode, transform::Bool)
    insupport(LengthDistribution(d), x) ? logpdf(d, x) : -Inf
end

function insupport(l::LengthDistribution, x::FNode)
    bl = get_branchlength_vector(x)
    all(isfinite.(bl)) && all(0.0 .< bl) && topo_placeholder(x, l) && !any(isnan.(bl))
end

#=
function insupport(d::BirthDeath, t::AbstractVector{T}) where {T<:Real}
    length(d) == length(t) && all(isfinite.(t)) && all(0 .< t)
end # function
=#

function insupport(t::UniformConstrained, x::FNode)::Bool
    topological.constraint_dict(t, x)
end

function insupport(t::UniformTopology, x::FNode)::Bool
    true
end

function topo_placeholder(x::FNode , l::LengthDistribution)
    true
end

function relistlength(d::CompoundDirichlet, x::AbstractArray)
  n = length(x)
  (Array(x), n)
end



"""
Using formula 11.24 with 11.25 from the following paper:
https://lukejharmon.github.io/pcm/chapter11_fitbd/#ref-FitzJohn2009-sg

"""
function _logpdf(d::BirthDeath, x::FNode)::Float64
    λ::Float64 = d.lambd
    μ::Float64 = d.mu
    f::Float64 = d.rho
    n::Int64 = length(get_leaves(x))
    h::Float64 = node_height(x)
    
    start::Float64 = λ ^ (n - 2)
    for node in post_order(x)
        node.root && continue
        #= 
        unsure about the t(k,t) & t(k,b) values in num & denum. The paper is not very clear 
        about those numbers, i.e. we don't know if they are both supposed to be positive, 
        negative or one positive & one negative
        =#
        num::Float64 = (f * λ - (μ - λ * (1 - f)) * exp((λ - μ) * 
                       node_height(node))) ^ 2
        denum::Float64 = (f * λ - (μ - λ * (1 - f)) * exp((λ - μ) *
                         (h - node_height(node)))) ^ 2
        value::Float64 = exp((λ - μ) * ((h - node_height(node)) - node_height(node))) *
                         (num / denum)
        start *= value
    end # for
    denumerator::Float64 = (1 - (1 - ((λ - μ) /(λ - (λ - μ) * exp((λ - μ) *
                           node_height(x)))))) ^ 2
    THEFORMULA = factorial(n - 1) * start / denumerator
end

"""
Fossilized Birth Death
Implemented following Heath, Huelsenbeck and Stadler 2013
DOI: 10.1073/pnas.1319091111
"""
mutable struct BirthDeathFossilized <: ContinuousMultivariateDistribution
    rho::Float64
    mu::Float64
    lambd::Float64
    psi::Float64
end # mutable struct

function _logpdf(d::FossilizedBirthDeath, x::FNode)::Float64

λ::Float64 = d.lambd
μ::Float64 = d.mu
ρ::Float64 = d.rho
ψ::Float64 = d.psi

c1::Float64 = sqrt((λ - μ - ψ)^2 + 4 * λ * ψ)
c2::Float64 = - ((λ - μ - 2 * λ *ρ - ψ )/ c1)
qt::Float64 = 2 * (1- c2 ^2) + exp(-c1) # I am not sure how to integrate t yet, need to read more of the paper
p0::Float64
p0_past::Float64
result::Float64 = 0.5 # placeholder

return result
end



"""
Strict Molecular Clock - Simplified Birth Death
Implemented folloing Yang & Rannala 1996
doi.org/10.1007/BF02338839
"""
mutable struct BirthDeathSimplified <: ContinuousMultivariateDistribution
    s::Int64
    mu::Float64
    lambd::Float64
end # mutable struct

length(d::BirthDeathSimplified) = d.s - 1

function insupport(d::BirthDeathSimplified, t::AbstractVector{T}) where {T<:Real}
    length(d) == length(t) && all(isfinite.(t)) && all(0 .< t)
end # function

function _logpdf(d::BirthDeathSimplified, t::AbstractVector{T}) where {T<:Real}

    f::Float64 = log(2.0)*(d.s-1.0)+log(d.mu)*(d.s-2.0)
    p0::Float = (d.s-2.0)*log((d.mu*(1.0-exp(-(d.lambd-d.mu)))) /
                (d.lambd-d.mu*exp(-(d.lambd-d.mu))))
    p0 += log(factorial(d.s)*(d.s-1.0))
    f -= po

    con::Float64 = 2.0*log(d.lambd - d.mu)

    for i in t
        f += (con-(d.lambd-d.mu)*i) - log(d.lambd-d.mu*exp(-(d.lambd-d.mu)*i))*2.0
    end # for
    return f
end # function
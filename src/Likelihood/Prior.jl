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


#################### Gradlogpdfs ####################

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

function gradlogpdf(d::BirthDeath, x::FNode)::Tuple{Float64, Vector{Float64}}
    blv::Vector{Float64} = get_branchlength_vector(x)[1:end]
    n::Int64 = length(get_leaves(x))
    p::Vector{FNode} = post_order(x)[1:end-1]
    sort!(p, by = x -> x.num)
    hv::Vector{Float64} = [n.height for n in p]
    g(y) = _logpdf(d, x, n, y, hv)
    r = Zygote.pullback(g, blv)
    r[1],r[2](1.0)[1]
end

function gradlogpdf(t::Union{UniformConstrained, UniformTopology, UniformBranchLength}, 
                    x::FNode)::Tuple{Float64, Vector{Float64}}

    blv = get_branchlength_vector(x)
    0.0, zeros(length(blv))
end

#################### Logpdfs ####################

function logpdf(d::CompoundDirichlet, x::FNode)
    internal_logpdf(d, get_branchlength_vector(x), internal_external(x))
end

function logpdf_sub(d::CompoundDirichlet, x::FNode, transform::Bool)
    insupport(LengthDistribution(d), x) ? logpdf(d, x) : -Inf
end

function logpdf(t::Union{UniformConstrained, UniformTopology, UniformBranchLength}, x::FNode)
    0.0
end

function logpdf(ex::exponentialBL, x::FNode)
    _logpdf(ex, x)
end

function _logpdf(d::exponentialBL, x::FNode)
    blv = get_branchlength_vector(x)
    sum(logpdf.(Exponential(d.scale), blv))
end

function logpdf(d::BirthDeath, x::FNode)
    blv::Vector{Float64} = get_branchlength_vector(x)
    n::Int64 = length(get_leaves(x))
    p::Vector{FNode} = post_order(x)[1:end-1]
    sort!(p, by = x -> x.num)
    hv::Vector{Float64} = [n.height for n in p]
    _logpdf(d, x, n, blv, hv)
end

"""
Using formula 11.24 with 11.25 from the following paper:
https://lukejharmon.github.io/pcm/chapter11_fitbd/#ref-FitzJohn2009-sg

The constant factorial factor (n-1)! is currently not in the formula, since the gradient
calculation with Zygote has problems with the BigInt number that it produces. Be aware of
"this multiplier if comparing likelihoods across different models for model selection" 
(from the paper above right before equation 11.18)

This paper & our method assume ultrametric trees.
"""
function _logpdf(d::BirthDeath, x::FNode, n::Int64, blv::Vector{Float64}, 
                 hv::Vector{Float64})::Float64
                
    λ::Float64 = d.lambd
    μ::Float64 = d.mu
    f::Float64 = d.rho
    h::Float64 = x.height
    
    start::Float64 = λ ^ (n - 2)

    for (i, bl) in enumerate(blv)
        num::Float64 = (f * λ - (μ - λ * (1 - f)) * exp((λ - μ) * hv[i])) ^ 2
        denum::Float64 = (f * λ - (μ - λ * (1 - f)) * exp((λ - μ) * (h - bl - hv[i]))) ^ 2
        value::Float64 = exp((λ - μ) * (h - bl - (2 * hv[i])) * (num / denum))
        # value::Float64 = exp((λ - μ) * (h - bl) * (num / denum))
        start *= value
    end # for

    denumerator::Float64 = (1 - (1 - ((λ - μ) /(λ - (λ - μ) * exp((λ - μ) * h))))) ^ 2
    log(start / denumerator)
end



#################### Insupports ####################

function insupport(l::LengthDistribution, x::FNode)
    bl = get_branchlength_vector(x)
    all(isfinite.(bl)) && all(0.0 .< bl) && topo_placeholder(x, l) && !any(isnan.(bl))
end

function insupport(d::BirthDeath, t::AbstractVector{T}) where {T<:Real}
    length(d) == length(t) && all(isfinite.(t)) && all(0 .< t)
end

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
Fossilized Birth Death
Implemented following Heath, Huelsenbeck and Stadler 2013
DOI: https://arxiv.org/pdf/1310.2968.pdf
"""
mutable struct BirthDeathFossilized <: ContinuousMultivariateDistribution
    rho::Float64
    lambd::Float64
    mu::Float64
    psi::Float64
    m::Int64
    k::Int64
end # mutable struct


function _logpdf(d::BirthDeathFossilized, x::FNode, internal_node_times::Vector{Float64}, 
                 leaf_fossil_times::Vector{Float64})::Float64

λ::Float64 = d.lambd
μ::Float64 = d.mu
ρ::Float64 = d.rho
ψ::Float64 = d.psi

start::Float64 = λ ^ (n + m - 1) * ψ ^ (k + m) * (4 * ρ) ^ n
c₁::Float64 = abs(sqrt((λ - μ - ψ)^2 + 4 * λ * ψ))
c₂::Float64 = - ((λ - μ - 2 * λ * ρ - ψ )/ c₁)

q(t) = 2 * (1 - c₂ ^ 2) + exp(-c₁ * t) * (1 - c₂) ^ 2  + exp(c₁ * t) * (1 + c₂) ^ 2
num(t)::Float64 = exp(-c₁ * t) * (1 - c₂) - (1 + c₂)
denum(t)::Float64 = exp(-c₁ * t) * (1 - c₂) + (1 + c₂) 
p₀(t)::Float64 = (λ + μ + ψ + c₁ * (num(t) / denum(t))) / (2 * λ)
p₁(t)::Float64 = (4 * ρ) / q(t)

for time in internal_node_times
    start *= 1 / q(time)
end

for time in leaf_fossil_times
    start *= p₀(time) * q(time)
 
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
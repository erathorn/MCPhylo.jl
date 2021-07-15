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

    r2 = first + second +third+fourth

    return r2

end

function gradlogpdf(d::CompoundDirichlet, x::T) where T <: GeneralNode
    int_ext = internal_external(x)
    blv = get_branchlength_vector(x)
    f(y) =  internal_logpdf(d, y, int_ext)
    r = Zygote.pullback(f, blv)
    return r[1],r[2](1.0)[1]
end

function gradlogpdf(t::Union{Missing, TopologyDistribution}, x::T
                   )::Tuple{Float64, Vector{Float64}} where T <: GeneralNode

    blv = get_branchlength_vector(x)
    0.0, zeros(length(blv))
end

function logpdf(d::CompoundDirichlet, x::T) where T <: GeneralNode
    internal_logpdf(d, get_branchlength_vector(x), internal_external(x))
end

function logpdf(t::Union{TopologyDistribution, Missing}, x::T) where T <: GeneralNode
    0.0
end

function logpdf(ex::exponentialBL, x::T) where T <: GeneralNode
    _logpdf(ex, x)
end

function insupport(l::LengthDistribution, x::T) where T <: GeneralNode
    bl = get_branchlength_vector(x)
    all(isfinite.(bl)) && all(0.0 .< bl) && topo_placeholder(x, l) && !any(isnan.(bl))
end 

function insupport(t::TopologyDistribution, x::T)::Bool where T <: GeneralNode
    topological(t.constraint_dict, x)
end

function logpdf_sub(d::CompoundDirichlet, x::T, transform::Bool) where T <: GeneralNode
    insupport(LengthDistribution(d), x) ? logpdf(d, x) : -Inf
end

function relistlength(d::CompoundDirichlet, x::AbstractArray)
  n = length(x)

  (Array(x), n)
end

function topo_placeholder(x::N , l::LengthDistribution) where N<:GeneralNode
    true
end

function _logpdf(d::exponentialBL, x::FNode)
    bl = get_branchlength_vector(x)
    sum(logpdf.(Exponential(d.scale), bl))
end

function gradlogpdf(d::exponentialBL, x::FNode)
    bl = get_branchlength_vector(x)
    g(y) = sum(logpdf.(Exponential(d.scale), y))
    r = Zygote.pullback(g, bl)
    r[1],r[2](1.0)[1]
end


"""
Strict Molecular Clock - BirthDeath
Implemented following Yang & Rannala 1997
doi.org/10.1093/oxfordjournals.molbev.a025811
"""
mutable struct BirthDeath <: ContinuousMultivariateDistribution
    s::Int64
    rho::Float64
    mu::Float64
    lambd::Float64
end # mutable struct

length(d::BirthDeath) = d.s - 1

function insupport(d::BirthDeath, t::AbstractVector{T}) where {T<:Real}
    length(d) == length(t) && all(isfinite.(t)) && all(0 .< t)
end # function

function _logpdf(d::BirthDeath, t::AbstractVector{T}) where {T<:Real}
    numerator::Float64 = (d.rho*(d.lambd-d.mu))/(d.rho*d.lambd + (d.lambd*(1.0-d.rho)-d.mu)*exp(d.mu-d.lambd))
    denum::Float64 = d.rho*(d.lambd-d.mu)
    vt1::Float64 = log(1.0-((denum*exp(d.mu-d.lambd))/(d.rho * numerator)))
    f::Float64 = log((2.0^(s-1.0))/(factorial(s)*(s-1.0)))
    for i in t
        f += ((d.lambd+(1.0-d.rho)+2.0*(denum-numerator))+(d.mu-d.lambd)*i)-vt1
    end # for
    return f
end # function

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
    p0::Float = (d.s-2.0)*log((d.mu*(1.0-exp(-(d.lambd-d.mu))))/(d.lambd-d.mu*exp(-(d.lambd-d.mu))))
    p0 += log(factorial(d.s)*(d.s-1.0))
    f -= po

    con::Float64 = 2.0*log(d.lambd - d.mu)

    for i in t
        f += (con-(d.lambd-d.mu)*i) - log(d.lambd-d.mu*exp(-(d.lambd-d.mu)*i))*2.0
    end # for
    return f
end # function

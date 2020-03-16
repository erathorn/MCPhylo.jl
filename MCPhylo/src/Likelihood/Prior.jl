
"""
    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64, nterm::Float64)
This structure implememts the CompoundDirichlet distribution described
in Zhang, Rannala and Yang 2012. (DOI:10.1093/sysbio/sys030)
"""
## Type declaration
mutable struct CompoundDirichlet <: ContinuousUnivariateDistribution
    alpha::Float64
    a::Float64
    beta::Float64
    c::Float64
    constraints::Union{Dict, Missing}

    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64) =
        new(alpha, a, beta, c, missing)

    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64, constraints::Dict) =
            new(alpha, a, beta, c, constraints)
end # struct

function internal_logpdf(d::CompoundDirichlet, b_lens::Any, int_leave_map::Vector{Int64}; rd::Bool=false)
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
    third = blen_leave_log*(d.a-1) + blen_int_log*(d.a*d.c-1.0)
    fourth = (d.alpha-d.a*nterm-d.a*d.c*n_int)*log(t_l)

    r2 = first + second +third+fourth

    return r2

end

function pgradient(d::CompoundDirichlet, x::Node)
    
    int_ext = internal_external(x)
    blv = get_branchlength_vector(x)


    g(mt) = internal_logpdf(d, mt, int_ext)

    r = val_der(g, blv)

    r[1], r[2][1]
end

function val_der(f, y...)
    gs1 = Flux.pullback(f, y...)
    gs1[1], gs1[2](1.0)
end


function _logpdf(d::CompoundDirichlet, x::Node)
    internal_logpdf(d, get_branchlength_vector(x), internal_external(x); rd=true)
end # function _logpdf

function insupport(d::CompoundDirichlet, x::Node)
    bl = get_branchlength_vector(x)
    all(isfinite.(bl)) && all(0.0 .<= bl) && topological(x, d.constraints) && !any(isnan.(bl))
end # function insupport


function logpdf_sub(d::ContinuousUnivariateDistribution, x::Node, transform::Bool)
    insupport(d, x) ? _logpdf(d, x) : -Inf
end


mutable struct exponentialBL <: ContinuousUnivariateDistribution
    scale::Float64
    constraints::Union{Dict, Missing}

    exponentialBL(scale::Float64) =
        new(scale, missing)
    exponentialBL(scale::Float64, c) =
            new(scale, c)
end

function _logpdf(d::exponentialBL, x::Node)
    bl = get_branchlength_vector(x)
    sum(bl)/d.scale - log(d.scale)*length(bl)
end

function insupport(d::exponentialBL, x::Node)
    bl = get_branchlength_vector(x)
    res = all(isfinite.(bl)) && all(0 .< bl) && topological(x, d.constraints) && !any(isnan.(bl))

    res
end # function insupport


function gradient(d::exponentialBL, x::Node)
    bl = get_branchlength_vector(x)
    ones(length(bl))./d.scale
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

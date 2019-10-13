
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

end # struct

length(d::CompoundDirichlet) = 23409
Base.size(d::CompoundDirichlet) = 1#(153,153)

function _logpdf(d::CompoundDirichlet, x::Node)

    xn = get_sum_seperate_length!(x)
    blen_int::Float64 = xn[1]
    blen_leave::Float64 = xn[2]
    t_l::Float64 = blen_int+blen_leave
    n_term::Float64 = length(get_leaves(x))
    n_int::Float64 = n_term-3.0
    ln1::Float64 = (d.a-1.0)*blen_leave + (d.a*d.c-1.0)*blen_int
    ln2::Float64 = (d.alpha-d.a*n_term-d.a*d.c*n_int)*log(t_l)- d.beta*t_l
    v1 = (d.alpha*log(d.beta))
    v2 = log(gamma(d.alpha))
    v3 = log(gamma(d.a*n_term+d.a*d.c*n_int))
    v4 = n_term*log(gamma(d.beta))
    v5 = n_int*log(gamma(d.beta*d.a))
    ln3::Float64 = (d.alpha*log(d.beta))-log(gamma(d.alpha))+log(gamma(d.a*n_term+d.a*d.c*n_int))-n_term*log(gamma(d.beta))-n_int*log(gamma(d.beta*d.a))
    return ln1+ln2+ln3
end # function _logpdf

function insupport(d::CompoundDirichlet, x::Node)

    all(isfinite.(get_branchlength_vector(x))) && all(0 .< get_branchlength_vector(x)) && topological(x, d.constraints)
end # function insupport


function logpdf_sub(d::ContinuousUnivariateDistribution, x::Node, transform::Bool)
    insupport(d, x) ? _logpdf(d, x) : -Inf
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

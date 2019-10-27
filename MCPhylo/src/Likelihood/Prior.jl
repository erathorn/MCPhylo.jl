
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

#length(d::CompoundDirichlet) = 23409
#Base.size(d::CompoundDirichlet) = 1#(153,153)

function internal_logpdf(d::CompoundDirichlet, b_lens::Any, int_leave_map::Array{Int64})
    blen_int = 0.0
    blen_leave = 0.0
    blen_int_log = 0.0
    blen_leave_log = 0.0
    nterm::Int64 = 0
    for i in eachindex(int_leave_map)
        if int_leave_map[i] == 1
            blen_int += b_lens[i]
            blen_int_log += log(b_lens[i])
        else
            blen_leave += b_lens[i]
            blen_leave_log += log(b_lens[i])
            nterm += 1
        end
    end

    t_l = blen_int+blen_leave
    n_int::Float64 = nterm-3.0

    first = (d.alpha*log(d.beta))-log(gamma(d.alpha)) - (t_l*d.beta)
    second = -log(gamma(d.a))-log(gamma(d.c))+log(gamma(d.a+d.c))
    third = blen_leave_log*(d.a-1) + blen_int_log*(d.a*d.c-1.0)
    fourth = (d.alpha-d.a*nterm-d.a*d.c*n_int)*log(t_l)

    r2 = first + second +third+fourth
    return r2

end

function gradient(d::CompoundDirichlet, x::Node)

    g(mt) = exp(internal_logpdf(d, mt, internal_external_map(x)))
    gradient(g, get_branchlength_vector(x))
end

function _logpdf(d::CompoundDirichlet, x::Node)

    internal_logpdf(d, get_branchlength_vector(x), internal_external_map(x))
    #xn = get_sum_seperate_length!(x)
    #blen_int::Float64 = xn[1]
    #blen_leave::Float64 = xn[2]
    #blen_int_log::Float64 = xn[3]
    #blen_leave_log::Float64 = xn[4]
    #t_l::Float64 = blen_int+blen_leave
    #n_term::Float64 = length(get_leaves(x))
    #n_int::Float64 = n_term-3.0

    #first = (d.alpha*log(d.beta))-log(gamma(d.alpha)) - (t_l*d.beta)
    #second = -log(gamma(d.a))-log(gamma(d.c))+log(gamma(d.a+d.c))
    #third = blen_leave_log*(d.a-1) + blen_int_log*(d.a*d.c-1.0)
    #fourth = (d.alpha-d.a*n_term-d.a*d.c*n_int)*log(t_l)

    #r2 = first + second +third+fourth


    #return r2
end # function _logpdf

function insupport(d::CompoundDirichlet, x::Node)
    bl = get_branchlength_vector(x)
    res = all(isfinite.(bl)) && all(0 .< bl) && topological(x, d.constraints) && !any(isnan.(bl))
    if !res
        println("insupport ",res)
        println(all(isfinite.(bl)) , all(0 .< bl), topological(x, d.constraints) ,!any(isnan.(bl)))
        println(bl)
    end
    res
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

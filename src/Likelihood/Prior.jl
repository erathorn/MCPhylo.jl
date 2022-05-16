function internal_logpdf(d::CompoundDirichlet, b_lens::Array{T},
                         int_leave_map::Vector{Int64})::T where T<:Real

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

function gradlogpdf(d::CompoundDirichlet, x::GeneralNode)
    int_ext = internal_external(x)
    blv = get_branchlength_vector(x)
    # use let block for proper capturing of variables
    f = let d=d, int_ext=int_ext
        y -> internal_logpdf(d, y, int_ext)
    end 
    r = Zygote.pullback(f, blv)
    return r[1],r[2](1.0)[1]
end

function gradlogpdf(d::exponentialBL, x::GeneralNode)
    bl = get_branchlength_vector(x)
        sum(logpdf.(Exponential(d.scale), bl)), -ones(length(bl))
end

function gradlogpdf(t::Union{UniformConstrained, UniformTopology, UniformBranchLength}, 
                    x::GeneralNode)::Tuple{Float64, Vector{Float64}}

    blv = get_branchlength_vector(x)
    0.0, zeros(length(blv))
end


function logpdf(d::CompoundDirichlet, x::GeneralNode)
    internal_logpdf(d, get_branchlength_vector(x), internal_external(x))
end

function logpdf(t::Union{UniformConstrained, UniformTopology, UniformBranchLength}, x::GeneralNode)
    0.0
end

function logpdf(ex::exponentialBL, x::GeneralNode)
    bl = get_branchlength_vector(x)
    sum(logpdf.(Exponential(ex.scale), bl))
end


function logpdf_sub(d::CompoundDirichlet, x::GeneralNode, transform::Bool)
    insupport(LengthDistribution(d), x) ? logpdf(d, x) : -Inf
end

function insupport(l::LengthDistribution, x::GeneralNode)
    bl = get_branchlength_vector(x)
    all(isfinite.(bl)) && topo_placeholder(x, l) && !any(isnan.(bl)) && all(0.0 .< bl)
     #&& 
end 

function insupport(t::UniformConstrained, x::GeneralNode)::Bool
    topological(t, x)
end

function insupport(t::UniformTopology, x::GeneralNode)::Bool
    true
end

function insupport(t::UniformTopology, x::AbstractArray)::Bool
    true
end

function topo_placeholder(x::GeneralNode , l::LengthDistribution)
    true
end

function relistlength(d::CompoundDirichlet, x::AbstractArray)
  n = length(x)
  (Array(x), n)
end


"""
Strict Molecular Clock - BirthDeath
Implemented following Yang & Rannala 1997
doi.org/10.1093/oxfordjournals.molbev.a025811
"""
struct BirthDeath <: ContinuousMultivariateDistribution
    s::Int64
    rho::Float64
    mu::Float64
    lambd::Float64
end # struct

length(d::BirthDeath) = d.s - 1

function insupport(d::BirthDeath, t::AbstractVector{T}) where {T<:Real}
    length(d) == length(t) && all(isfinite.(t)) && all(0 .< t)
end # function

function _logpdf(d::BirthDeath, t::AbstractVector{T}) where {T<:Real}
    numerator::Float64 = (d.rho*(d.lambd-d.mu))/(d.rho*d.lambd + 
                         (d.lambd*(1.0-d.rho)-d.mu)*exp(d.mu-d.lambd))
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
struct BirthDeathSimplified <: ContinuousMultivariateDistribution
    s::Int64
    mu::Float64
    lambd::Float64
end # struct

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
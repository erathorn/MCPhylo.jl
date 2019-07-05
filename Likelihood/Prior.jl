module Prior

using Markdown
using SpecialFunctions
using Distributions


import Distributions: length, insupport, _logpdf

"""
    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64, nterm::Float64)
This structure implememts the CompoundDirichlet distribution described
in Zhang, Rannala and Yang 2012. (DOI:10.1093/sysbio/sys030)
"""
## Type declaration
mutable struct CompoundDirichlet <: ContinuousMultivariateDistribution
    alpha::Float64
    a::Float64
    beta::Float64
    c::Float64
    n_term::Float64
end # struct

length(d::CompoundDirichlet) = 2# * d.n_term - 3.0

function _logpdf(d::CompoundDirichlet, x::AbstractVector{T}) where {T<:Real}
    blen_int::Float64 = x[1]
    blen_leave::Float64 = x[2]
    t_l::Float64 = blen_int+blen_leave
    n_int::Float64 = d.n_term-3.0
    ln1::Float64 = (d.a-1.0)*blen_leave + (d.a*d.c-1.0)*blen_int
    ln2::Float64 = (d.alpha-d.a*d.n_term-d.a*d.c*n_int)*log(t_l)- d.beta*t_l
    ln3::Float64 = (d.alpha*log(d.beta))-log(gamma(d.alpha))+log(gamma(d.a*d.n_term+d.a*d.c*n_int))-d.n_term*log(gamma(d.beta))-n_int*log(gamma(d.beta*d.a))
    return ln1+ln2+ln3
end # function _logpdf

function insupport(d::CompoundDirichlet, x::AbstractVector{T}) where {T<:Real}
    length(d) == length(x) && all(isfinite.(x)) && all(0 .<= x)
end # function insupport

end  # module Prior


mutable struct MultivariateUniformTrunc <: ContinuousMultivariateDistribution
  dims::Int
  lower::Float64
  upper::Float64
end

Base.size(d::MultivariateUniformTrunc) = (d.dims, )
length(d::MultivariateUniformTrunc) = d.dims

function insupport(d::MultivariateUniformTrunc, x::Vector{Float64})

    all(isfinite.(x[1:end-1])) && all(d.lower .< x[1:end-1] .< d.upper)
end # function insupport

function rand(d::MultivariateUniformTrunc)
  [rand(Uniform(d.lower, d.upper)) for i in 1:d.dims]
end

function logpdf(d::MultivariateUniformTrunc, x::Vector{Float64})
  o = 0.0
  for i in eachindex(x[1:end-1])
    o += logpdf(Uniform(d.lower, d.upper), x[i])
  end
  o
end



mutable struct CompoundDirichletWrap <: ContinuousMultivariateDistribution
  dims::Int
  alpha::Float64
  a::Float64
  beta::Float64
  c::Float64
  map::Vector{Int64}

  CompoundDirichletWrap(n::N, alpha::Float64,a::Float64,beta::Float64,c::Float64) where {N<:TreeVariate} =
    new(length(get_branchlength_vector(n)), alpha,a, beta, c, internal_external_map(n))
end


Base.size(d::CompoundDirichletWrap) = (d.dims, )
length(d::CompoundDirichletWrap) = d.dims

function insupport(d::CompoundDirichletWrap, x::Vector{Float64})
    all(isfinite.(x[1:end-1])) && all(0 .< x[1:end-1])
end # function insupport


function logpdf(d::CompoundDirichletWrap, x::Vector{Float64})
  blen_int::Float64 = 0
  blen_leave::Float64 = 0
  n_term::Float64 = 0
  for i in eachindex(d.map)
    if d.map[i] == 0
      blen_leave += x[i]
      n_term += 1.0
    else
      blen_int += x[i]
    end
  end

    t_l::Float64 = blen_int+blen_leave
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
end

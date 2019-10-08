
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

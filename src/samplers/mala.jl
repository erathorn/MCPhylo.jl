#################### Metropolis-Adjusted Langevin Algorithm ####################

#################### Types ####################

mutable struct MALATune <: SamplerTune
  logfgrad::Union{Function, Missing}
  epsilon::Float64
  SigmaL::Union{UniformScaling{Bool}, LowerTriangular{Float64}}

  MALATune() = new()

  MALATune(x::Vector, epsilon::Real) = new(missing, epsilon, I)

  MALATune(x::Vector, epsilon::Real, logfgrad::Function) =
    new(logfgrad, epsilon, I)

  MALATune(x::Vector, epsilon::Real, Sigma::Matrix{T}) where {T<:Real} =
    new(missing, epsilon, cholesky(Sigma).L)

  function MALATune(x::Vector, epsilon::Real, Sigma::AbstractMatrix{T},
                    logfgrad::Function) where {T<:Real}
    new(logfgrad, epsilon, cholesky(Sigma).L)
  end
end


const MALAVariate = SamplerVariate{MALATune}

validate(v::MALAVariate) = validate(v, v.tune.SigmaL)

validate(v::MALAVariate, SigmaL::UniformScaling) = v

function validate(v::MALAVariate, SigmaL::LowerTriangular)
  n = length(v)
  size(SigmaL, 1) == n ||
    throw(ArgumentError("Sigma dimension differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################
"""
    MALA(params::ElementOrVector{Symbol}, epsilon::Real; args...)

Construct a `Sampler` object for MALA sampling. Parameters are assumed to be
continuous, but may be constrained or unconstrained.

Returns a `Sampler{MALATune}`` type object.

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `epsilon`: factor by which the drift and covariance matrix of the proposal distribution are scaled.

* `Sigma`: covariance matrix for the multivariate normal proposal distribution. The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated. If omitted, the identity matrix is assumed.

* `dtype` : differentiation for gradient calculations. Options are
    * `:central` : central differencing
    * `:forward` : forward differencing.
"""
function MALA(params::ElementOrVector{Symbol}, epsilon::Real; args...)
  MALASampler(params, epsilon; args...)
end
"""
    MALA(params::ElementOrVector{Symbol}, epsilon::Real, Sigma::Matrix{T}; args...)

Construct a `Sampler` object for MALA sampling. Parameters are assumed to be
continuous, but may be constrained or unconstrained.

Returns a `Sampler{MALATune}`` type object.
"""
function MALA(params::ElementOrVector{Symbol}, epsilon::Real,
               Sigma::Matrix{T}; args...) where {T<:Real}
  MALASampler(params, epsilon, Sigma; args...)
end

function MALASampler(params, pargs...; dtype::Symbol=:forward)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, pargs...)
    sample!(v, x -> logpdfgrad!(block, x, dtype))
    relist(block, v)
  end
  Sampler(params, samplerfx, MALATune())
end


#################### Sampling Functions ####################

sample!(v::MALAVariate) = sample!(v, v.tune.logfgrad)
"""
    sample!(v::MALAVariate, logfgrad::Function)

Draw one sample from a target distribution using the MALA sampler. Parameters
are assumed to be continuous and unconstrained.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::MALAVariate, logfgrad::Function)
  tune = v.tune

  L = sqrt(tune.epsilon) * tune.SigmaL
  Linv = inv(L)
  M2 = 0.5 * L * L'

  logf0, grad0 = logfgrad(v.value)
  y = v + M2 * grad0 + L * randn(length(v))
  logf1, grad1 = logfgrad(y)

  q0 = -0.5 * sum(abs2, Linv * (v - y - M2 * grad1))
  q1 = -0.5 * sum(abs2, Linv * (y - v - M2 * grad0))

  if rand() < exp((logf1 - q1) - (logf0 - q0))
    v[:] = y
  end

  v
end

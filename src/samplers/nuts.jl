#################### No-U-Turn Sampler ####################

#################### Types and Constructors ####################

mutable struct NUTSTune <: SamplerTune
  logf::Union{Function, Missing}
  adapt::Bool
  alpha::Float64
  epsilon::Float64
  epsilonbar::Float64
  gamma::Float64
  Hbar::Float64
  kappa::Float64
  m::Int
  mu::Float64
  nalpha::Int
  t0::Float64
  target::Float64
  tree_depth::Int

  NUTSTune() = new()

  function NUTSTune(x::Vector, epsilon::Real, logfgrad::Union{Function, Missing};
                    target::Real=0.6, tree_depth::Int=10)
    new(logfgrad, false, 0.0, epsilon, 1.0, 0.05, 0.0, 0.75, 0, NaN, 0, 10.0,
        target, tree_depth)
  end
end

#NUTSTune(x::Vector, epsilon::Real; args...) =
#  NUTSTune(x, epsilon; args...)


const NUTSVariate = Sampler{NUTSTune, T} where T


#################### Sampler Constructor ####################
"""
    NUTS(params::ElementOrVector{Symbol}; dtype::Symbol=:forward, args...)

Construct a `Sampler` object for NUTS sampling, with the algorithm’s step size
parameter adaptively tuned during burn-in iterations. Parameters are assumed
to be continuous, but may be constrained or unconstrained.

Returns a `Sampler{NUTSTune}` type object.

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `dtype` : differentiation for gradient calculations. Options are
    * `:central` : central differencing
    * `:forward` : forward differencing.

* `args...`: additional keyword arguments to be passed to the `NUTSVariate` constructor.
"""
function NUTS(params::ElementOrVector{Symbol}; epsilon::Real = -Inf, kwargs...)
  tune = NUTSTune(Float64[], epsilon, logpdfgrad!; kwargs...)
  Sampler(Float64[], params, tune, Symbol[], true)
  #Sampler(params, samplerfx, NUTSTune())
end


#################### Sampling Functions ####################

sample!(v::NUTSVariate; args...) = sample!(v, v.tune.logfgrad; args...)
"""
    sample!(v::NUTSVariate, logfgrad::Function; adapt::Bool=false)

Draw one sample from a target distribution using the NUTS sampler. Parameters
are assumed to be continuous and unconstrained.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::NUTSVariate{T}, logfgrad::Function; adapt::Bool=false) where T<: AbstractArray{<: Real}
  tune = v.tune
  
  if tune.m == 0 && isinf(tune.epsilon)
    tune.epsilon = nutsepsilon(v.value, logfgrad, tune.target)
  end
  setadapt!(v, adapt)
  if tune.adapt
    tune.m += 1
    nuts_sub!(v, tune.epsilon, logfgrad)
    p = 1.0 / (tune.m + tune.t0)
    tune.Hbar = (1.0 - p) * tune.Hbar +
                p * (tune.target - tune.alpha / tune.nalpha)
    tune.epsilon = exp(tune.mu - sqrt(tune.m) * tune.Hbar / tune.gamma)
    p = tune.m^-tune.kappa
    tune.epsilonbar = exp(p * log(tune.epsilon) +
                          (1.0 - p) * log(tune.epsilonbar))
  else
    if (tune.m > 0) tune.epsilon = tune.epsilonbar end
    nuts_sub!(v, tune.epsilon, logfgrad)
  end
  v
end


function setadapt!(v::NUTSVariate{T}, adapt::Bool) where T<: AbstractArray{<: Real}
  tune = v.tune
  if adapt && !tune.adapt
    tune.m = 0
    tune.mu = log(10.0 * tune.epsilon)
  end
  tune.adapt = adapt
  v
end


function nuts_sub!(v::NUTSVariate{T}, epsilon::Real, logfgrad::Function) where T<: AbstractArray{<: Real}
  n = length(v)
  x = deepcopy(v.value)
  logf, grad = logfgrad(x)
  r = randn(n)
  #x, r, logf, grad = leapfrog(v.value, randn(n), zeros(n), 0.0, logfgrad)
  logp0 = logf - 0.5 * dot(r, r)
  logu0 = logp0 + log(rand())
  xminus = xplus = deepcopy(x)
  rminus = rplus = deepcopy(r)
  gradminus = gradplus = deepcopy(grad)
  j = 0
  n = 1
  s = true
  while s && j < v.tune.tree_depth
    pm = 2 * (rand() > 0.5) - 1
    if pm == -1
      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha,
        nalpha = buildtree(xminus, rminus, gradminus, pm, j, epsilon, logfgrad,
                           logp0, logu0)
    else
      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha =
        buildtree(xplus, rplus, gradplus, pm, j, epsilon, logfgrad, logp0,
                  logu0)
    end
    if sprime && rand() < nprime / n
      v[:] = xprime
    end
    j += 1
    n += nprime
    s = sprime && nouturn(xminus, xplus, rminus, rplus)
    v.tune.alpha, v.tune.nalpha = alpha, nalpha
  end
  v
end


function leapfrog(x::Vector{Float64}, r::Vector{Float64}, grad::Vector{Float64},
                  epsilon::Real, logfgrad::Function)
  r += (0.5 * epsilon) * grad
  x += epsilon * r
  logf, grad = logfgrad(x)
  r += (0.5 * epsilon) * grad
  x, r, logf, grad
end


function buildtree(x::Vector{Float64}, r::Vector{Float64},
                   grad::Vector{Float64}, pm::Integer, j::Integer,
                   epsilon::Real, logfgrad::Function, logp0::Real, logu0::Real)
  if j == 0
    xprime, rprime, logfprime, gradprime = leapfrog(x, r, grad, pm * epsilon,
                                                    logfgrad)
    logpprime = logfprime - 0.5 * dot(rprime, rprime)
    nprime = Int(logu0 < logpprime)
    sprime = logu0 < logpprime + 1000.0
    xminus = xplus = xprime
    rminus = rplus = rprime
    gradminus = gradplus = gradprime
    alphaprime = min(1.0, exp(logpprime - logp0))
    nalphaprime = 1
  else
    xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
      alphaprime, nalphaprime = buildtree(x, r, grad, pm, j - 1, epsilon,
                                          logfgrad, logp0, logu0)
    if sprime
      if pm == -1
        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 = buildtree(xminus, rminus, gradminus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 = buildtree(xplus, rplus, gradplus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0)
      end
      if rand() < nprime2 / (nprime + nprime2)
        xprime = xprime2
      end
      nprime += nprime2
      sprime = sprime2 && nouturn(xminus, xplus, rminus, rplus)
      alphaprime += alphaprime2
      nalphaprime += nalphaprime2
    end
  end
  xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
    alphaprime, nalphaprime
end


function nouturn(xminus::Vector{Float64}, xplus::Vector{Float64},
                 rminus::Vector{Float64}, rplus::Vector{Float64})
  xdiff = xplus - xminus
  dot(xdiff, rminus) >= 0 && dot(xdiff, rplus) >= 0
end

function nutsepsilon(x::Vector{Float64}, logfgrad::Function, target::Float64)
  n = length(x)
  _, r0, logf0, grad0 = leapfrog(x, randn(n), zeros(n), 0.0, logfgrad)
  epsilon = 1.0
  _, rprime, logfprime, gradprime = leapfrog(x, r0, grad0, epsilon, logfgrad)
  prob = exp(logfprime - logf0 - 0.5 * (dot(rprime, rprime) - dot(r0, r0)))
  pm = 2 * (prob > target) - 1
  while prob^pm > target^pm
    epsilon *= 2.0^pm
    _, rprime, logfprime, _ = leapfrog(x, r0, grad0, epsilon, logfgrad)
    prob = exp(logfprime - logf0 - 0.5 * (dot(rprime, rprime) - dot(r0, r0)))
  end
  epsilon
end
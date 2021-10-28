#################### Binary Metropolised Gibbs Sampler ####################

#################### Types and Constructors ####################

const BMGForm = Union{Int, Vector{Vector{Int}}}

mutable struct BMGTune{F<:BMGForm} <: SamplerTune
  logf::Union{Function, Missing}
  k::F

  BMGTune{F}() where {F<:BMGForm} = new{F}()

  BMGTune{F}(x::Vector, k::F) where {F<:BMGForm} = new{F}(missing, k)

  BMGTune{F}(x::Vector, k::F, logf::Function) where {F<:BMGForm} =
    new{F}(logf, k)
end


const BMGIntVariate = Sampler{BMGTune{Int}, T} where T
const BMGVecVariate = Sampler{BMGTune{Vector{Vector{Int}}}, T} where T

BMGVariate(x::Vector, logf::Function; k::F=1) where {F<:BMGForm} =
  Sampler{BMGTune{F}}(x, k, logf)


function validate(v::BMGIntVariate)
  n = length(v)
  v.tune.k <= n || throw(ArgumentError("k exceeds variate length $n"))
  validatebinary(v)
end

function validate(v::BMGVecVariate)
  n = length(v)
  mapreduce(maximum, max, v.tune.k) <= n ||
    throw(ArgumentError("indices in k exceed variate length $n"))
  validatebinary(v)
end


#################### Sampler Constructor ####################
"""
    BMG(params::ElementOrVector{Symbol}; k::F=1) where {F<:BMGForm}

Construct a `Sampler` object for BMG sampling. Parameters are assumed to have
 binary numerical values (0 or 1).

Returns a `Sampler{BMGTune{typeof(k)}}` type object.

* `params`: stochastic node(s) to be updated with the sampler.

* `k`:  number of parameters or vector of parameter indices to select at random for simultaneous updating in each call of the sampler.
"""
function BMG(params::ElementOrVector{Symbol}; k::F=1) where {F<:BMGForm}
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = Sampler(block, k)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMGTune{F}())
end


#################### Sampling Functions ####################

sample!(v::Sampler{BMGTune{F}}) where {F<:BMGForm} = sample!(v, v.tune.logf)
"""
    sample!(v::Sampler{BMGTune{F}}, logf::Function) where {F<:BMGForm}

Draw one sample from a target distribution using the BMG sampler. Parameters
are assumed to have binary numerical values (0 or 1).

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::Sampler{BMGTune{F}}, logf::Function) where {F<:BMGForm}
  n = length(v)
  probs = Vector{Float64}(undef, n)
  idx = randind(v)

  pbernoulli! = function(x, probs)
    for i in idx
      x_i = x[i]

      x[i] = 0.0
      logf0 = logf(x)
      x[i] = 1.0
      logf1 = logf(x)

      p = invlogit(logf1 - logf0)
      probs[i] = 0.0 < p < 1.0 ? p : 0.5

      x[i] = x_i
    end
    probs
  end

  x = v[:]
  pbernoulli!(x, probs)
  theta = map(p -> rand() < p, probs[idx])

  if n == 1
    v[idx] = theta
  else
    y = v[:]
    y[idx] = theta
    qy = mapreduce(i -> y[i] == 1.0 ? log(probs[i]) : log(1.0 - probs[i]), +,
                   idx)

    pbernoulli!(y, probs)
    qx = mapreduce(i -> x[i] == 1.0 ? log(probs[i]) : log(1.0 - probs[i]), +,
                   idx)

    if rand() < exp((logf(y) - qy) - (logf(x) - qx))
      v[idx] = theta
    end
  end

  v
end

randind(v::BMGIntVariate) = sample(1:length(v), v.tune.k, replace=false)
randind(v::BMGVecVariate) = sample(v.tune.k)

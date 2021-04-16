#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

mutable struct RWMCTune <: SamplerTune
  logf::Union{Function, Missing}

  RWMCTune() = new()

  #function RWMCTune(x::Vector, logf::Union{Function, Missing}
  #                 )
  #  new(logf)
  #end

  function RWMCTune(x::Vector,logf::Union{Function, Missing})
    new(logf)
  end
end

RWMCTune(x::Vector; args...) =
  RWMCTune(x, missing; args...)


const RWMCVariate = SamplerVariate{RWMCTune}

validate(v::RWMCVariate) = v#validate(v, v.tune.scale)

#validate(v::RWMCVariate, scale::Float64) = v

#function validate(v::RWMCVariate, scale::Vector)
#  n = length(v)
#  length(scale) == n ||
#    throw(ArgumentError("length(scale) differs from variate length $n"))
#  v
#end


#################### Sampler Constructor ####################

function RWMC(params::ElementOrVector{Symbol};
             args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, RWMCTune())
end


#################### Sampling Functions ####################

sample!(v::RWMCVariate) = sample!(v, v.tune.logf)

function transformfunc(x, in_min, in_max, out_min, out_max)
  return trunc(Int64, (x-in_min)*(out_max-out_min) / (in_max-in_min) + out_min)
end

function sample!(v::RWMCVariate, logf::Function)

  x = v + map(x->transformfunc(x, 1, 3, -1, 1), rand(Categorical(ones(3)./3), length(v)))
  #x = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end

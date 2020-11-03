#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################

mutable struct DMHTune <: SamplerTune
  outer<:SamplerVariate
  inner<:SamplerVariate

  DMHTune() = new()
  function DMHTune(outer::S, inner::T) where {S<:SamplerVariate, T<:SamplerVariate}
    typeof(outer) == DMHVariate && throw(ArgumentError($outer " cannot be of type " $typeof(outer)))
    typeof(inner) == DMHVariate && throw(ArgumentError($inner " cannot be of type " $typeof(inner)))
    new(outer, inner)
  end
end

const DMHVariate = SamplerVariate{DMHTune}


#################### Sampler Constructor ####################

function DMH(params::ElementOrVector{Symbol}, outer::Symbol,
             inner::Symbol)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    f = (x, sz, ll, gr) -> DMH_logpdf!(block, x, sz, ll, gr)
    v = SamplerVariate(block, f, NullFunction(); args...)

    DMH_sample!(v::DMHVariate, f)

    relist(block, v)
  end
  Sampler(params, samplerfx, DMHTune())
end

function DMH_sample!(v::DMHVariate, args...)
  # 1. propose theta prime
  # 2. Generate the auxillary variable using theta prime

  # this way of generating theta_prime from the current values of theta
  # taes care of teh transition probability from theta_prim to theta and vice versa
  # the values equal and will cancel out.
  x = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))



end

#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

mutable struct RWMTune <: SamplerTune
  logf::Union{Function, Missing}
  scale::Union{Float64, Vector{Float64}}
  proposal::SymDistributionType

  RWMTune() = new()

  function RWMTune(x::Vector, scale::Real, logf::Union{Function, Missing};
                   proposal::SymDistributionType=Normal)
    new(logf, Float64(scale), proposal)
  end

  function RWMTune(x::Vector, scale::Vector{T},
                  logf::Union{Function, Missing};
                  proposal::SymDistributionType=Normal) where {T<:Real}
    new(logf, convert(Vector{Float64}, scale), proposal)
  end
end

RWMTune(x::Vector, scale::ElementOrVector{T}; args...) where {T<:Real} =
  RWMTune(x, scale, missing; args...)


const RWMVariate = SamplerVariate{RWMTune}

validate(v::RWMVariate) = validate(v, v.tune.scale)

validate(v::RWMVariate, scale::Float64) = v

function validate(v::RWMVariate, scale::Vector)
  n = length(v)
  length(scale) == n ||
    throw(ArgumentError("length(scale) differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################
"""
    RWM(params::ElementOrVector{Symbol},
                  scale::ElementOrVector{T}; args...) where {T<:Real})

Construct a `Sampler` object for RWM sampling. Parameters are assumed to be
continuous, but may be constrained or unconstrained.

Returns a `Sampler{RWMTune}` type object.
"""
function RWM(params::ElementOrVector{Symbol},
              scale::ElementOrVector{T}; args...) where {T<:Real}

  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, scale; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, RWMTune())
end

"""
  RWM(params::ElementOrVector{Symbol}, moves::Array{Symbol}; args...)

Construct the RWM sampler for Trees. If you set moves to :all it will use all
eligible moves to change the tree. These are currently:
NNI, Slide, Swing, :EdgeLength
"""
function RWM(params::ElementOrVector{Symbol}, moves::ElementOrVector{Symbol}; args...) where {T<:Real}
  #eligible = [:NNI, :SPR, :Slide, :Swing, :EdgeLength] # use after SPR branch is merged
  eligible = [:NNI, :Slide, :Swing, :EdgeLength]
  to_use = Symbol[]
  if moves == :all
    to_use = eligible
  else
    for i in moves
      !(i in eligible) &&
        throw("$i is not an eligible tree move. The list of eligible tree moves is $eligible")
      push!(to_use, i)
    end
  end
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    # scale is not necessary for RWM on trees so just set it to 1
    v = SamplerVariate(block, 1.0; args...)
    sample!(v, x -> logpdf!(block, x), to_use)
    relist(block, v)
  end
  Sampler(params, samplerfx, RWMTune())
end



#################### Sampling Functions ####################

sample!(v::RWMVariate) = sample!(v, v.tune.logf)



function sample!(v::RWMVariate, logf::Function, moves::Array{Symbol})
  tree = v[1]
  tc = deepcopy(tree)
  move = rand(moves)
  if move == :NNI
    NNI!(tree)
  elseif move == :SPR
    tree = SPR(tree)
  elseif move == :Slide
    slide!(tree)
  elseif move == :Swing
    swing!(tree)
  elseif move == :EdgeLength
    change_edge_length!(tree)
  else
    throw("Tree move not elegible ")
  end
  if rand() < exp(logf(tree) - logf(tc))
    v[1] = tree
  else
    v[1] = tc
  end
  v
end

"""
    sample!(v::RWMVariate, logf::Function)

Draw one sample from a target distribution using the RWM sampler. Parameters
are assumed to be continuous and unconstrained.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::RWMVariate, logf::Function)
  x = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end

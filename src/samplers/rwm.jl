#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

mutable struct RWMTune <: SamplerTune
  logf::Union{Function, Missing}
  scale::Union{Float64, Vector{Float64}}
  eligible::Vector{Symbol}
  proposal::SymDistributionType

  RWMTune() = new()

  function RWMTune(x::Vector, scale::Real, logf::Union{Function, Missing}, eligible::Vector{Symbol};
                   proposal::SymDistributionType=Normal)
    new(logf, Float64(scale), eligible, proposal)
  end

  function RWMTune(x::Vector, scale::Vector{T},
                  logf::Union{Function, Missing},
                  eligible::Vector{Symbol};
                  proposal::SymDistributionType=Normal) where {T<:Real}
    new(logf, convert(Vector{Float64}, scale), eligible, proposal)
  end
end

const RWMVariate = SamplerVariate{RWMTune, T} where T

validate(v::RWMVariate{T}) where T<:AbstractArray{S} where S<:Union{Real, GeneralNode} = validate(v, v.tune.scale)

validate(v::RWMVariate{T}, scale::Float64) where T<:AbstractArray{S} where S<:Union{Real, GeneralNode} = v 

function validate(v::RWMVariate{T}, scale::Vector) where T<:AbstractArray{S} where S<:Union{Real, GeneralNode}
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

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `scale`: scaling value or vector of the same length as the combined elements of nodes `params` for the `proposal` distribution. Values are relative to the unconstrained parameter space, where candidate draws are generated.

* `args...`: additional keyword arguments to be passed to the `RWMVariate` constructor.
"""
function RWM(params::ElementOrVector{Symbol},
              scale::ElementOrVector{T}; args...) where {T<:Real}
  tune = RWMTune(Float64[], scale, logpdf!, Symbol[])
  Sampler(params, tune, Symbol[], false)
end

"""
    RWM(params::ElementOrVector{Symbol}, moves::Array{Symbol}; args...)

Construct the RWM sampler for Trees. If you set moves to :all it will use all
eligible moves to change the tree. These are currently:
NNI, SPR, Slide, Swing, :EdgeLength

Returns a `Sampler{RWMTune}` type object.
"""
function RWM(params::ElementOrVector{Symbol}, moves::ElementOrVector{Symbol}; args...)
  eligible = [:NNI, :SPR, :Slide, :Swing, :EdgeLength]
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
  
  tune = RWMTune(Float64[], logpdf!, 1.0, to_use)
  Sampler(params, tune, Symbol[], transform)
end



#################### Sampling Functions ####################

sample!(v::RWMVariate) = sample!(v, v.tune.logf)


"""
    sample!(v::RWMVariate, logf::Function, moves::Array{Symbol})

Propose a new tree by randomly performing a move from the ones specified in `moves`.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::RWMVariate{T}, logf::Function; kwargs...) where T<:GeneralNode
  tree = v[1]
  tc = deepcopy(tree)

  move = rand(v.elegible)
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
function sample!(v::RWMVariate{T}, logf::Function; kwargs...) where T<:AbstractArray{<:Real}
  x = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end

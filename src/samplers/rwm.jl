#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

mutable struct RWMTune{F<:Function,T<:Union{Float64,Vector{Float64}}, D<:SymDistributionType} <: SamplerTune
    logf::F
    scale::T
    eligible::Vector{Symbol}
    proposal::D

    function RWMTune(
        x::Vector,
        scale::Real,
        logf::F,
        eligible::Vector{Symbol};
        proposal::D = Normal,
    ) where {F, D}
        new{F,Float64, Type{proposal}}(logf, Float64(scale), eligible, proposal)
    end

    function RWMTune(
        x::Vector,
        scale::Vector{T},
        logf::F,
        eligible::Vector{Symbol};
        proposal::D = Normal,
    ) where {T<:Real,F, D}
        new{F,Vector{Float64}, Type{proposal}}(logf, convert(Vector{Float64}, scale), eligible, proposal)
    end
end

const RWMVariate = Sampler{RWMTune{F,S,D},T} where {T,F,S,D}

validate(v::Sampler{RWMTune{F,S,D},T}) where {F,S,T,D} = validate(v, v.tune.scale)

validate(v::Sampler{RWMTune{F,S,D},T}, scale::Float64) where {F,S,T,D} = v

function validate(v::Sampler{RWMTune{F,S,D},T}, scale::Vector) where {F,S,T,D}
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
function RWM(
    params::ElementOrVector{Symbol},
    scale::ElementOrVector{T};
    proposal = Normal,
    args...,
) where {T<:Real}
    prop_dist =  proposal
    tune = RWMTune(Float64[], scale, logpdf!, Symbol[], proposal=prop_dist)
    Sampler(params, tune, Symbol[], true)
end

"""
    RWM(params::ElementOrVector{Symbol}, moves::Array{Symbol}; args...)

Construct the RWM sampler for Trees. If you set moves to :all it will use all
eligible moves to change the tree. These are currently:
NNI, SPR, Slide, Swing, :EdgeLength

Returns a `Sampler{RWMTune}` type object.
"""
function RWM(params::ElementOrVector{Symbol}, moves::ElementOrVector{Symbol}; kwargs...)
    eligible = [:NNI, :SPR, :Slide, :Swing, :EdgeLength]
    to_use = Symbol[]
    if moves == :all
        to_use = eligible
    else
        for i in moves
            !(i in eligible) && throw(
                "$i is not an eligible tree move. The list of eligible tree moves is $eligible",
            )
            push!(to_use, i)
        end
    end

    tune = RWMTune(Float64[], 1.0, logpdf!, to_use)
    Sampler(params, tune, Symbol[], false)
end



#################### Sampling Functions ####################


"""
    sample!(v::RWMVariate, logf::Function, moves::Array{Symbol})

Propose a new tree by randomly performing a move from the ones specified in `moves`.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(
    v::Sampler{RWMTune{F,S,D},Vector{T}},
    logf::Function;
    kwargs...,
) where {T<:GeneralNode,F,S,D}
    tree = v[1]
    tc = deepcopy(tree)

    move = rand(v.tune.eligible)
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
    if rand() < exp(logf([tree]) - logf([tc]))
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
function sample!(
    v::Sampler{RWMTune{F,S,D},T},
    logf::Function;
    kwargs...,
) where {T<:AbstractArray{<:Real},F,S,D}
    x = v + propose(v.tune.proposal(), length(v), v.tune.scale)
    if rand() < exp(logf(x) - logf(v.value))
        v[:] = x
    end
    v
end

function propose(probDist::D, n::Int, scale::ElementOrVector{Float64}) where D<:ContinuousDistribution
    scale .* rand(probDist, n)
end

function propose(probDist::D, n::Int, scale::ElementOrVector{Float64}) where D<:DiscreteDistribution
    rand(probDist, n)
end

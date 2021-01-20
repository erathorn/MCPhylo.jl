#################### Empirical Move ####################

#################### Types and Constructors ####################


mutable struct EmpiricalTune <: SamplerTune
    samplerfun::Union{Function, Missing}
    width::Int64
    replacement::Bool

    EmpiricalTune() = new()

    EmpiricalTune(samplerfun::Union{Function, Missing}, width::Int64, replacement::Bool=false)=
        new(samplerfun, width, false)
end

EmpiricalTune(x::Vector, width::Int64, replacement::Bool=false)=
    EmpiricalTune(missing, width, replacement)

const EmpiricalVariate = SamplerVariate{EmpiricalTune}


#################### Sampler Constructor ####################
"""

function Empirical(params::ElementOrVector{Symbol}, width::Int64; args...)

This Sampler generates Samples from a distribution by selecting a random value
from this distribution. The move will always be accepted.
"""
function Empirical(params::ElementOrVector{Symbol},
                width::Int64;
                args...)
    samplerfx = function(model::Model, block::Integer)
        block = SamplingBlock(model, block)
        v = SamplerVariate(block, width, args...)
        sample!(v, x -> rand!(block, x))
        relist(block, v)
    end
    Sampler(params, samplerfx, EmpiricalTune())
end

sample!(v::EmpiricalVariate) = sample(v, v.tune.samplerfun)
function sample!(v::EmpiricalVariate, sampler_func::Function)
    width = v.tune.width
    sam = sampler_func(width)
    if width > 1 && v.tune.replacement
        while length(Set(sam)) != width
            sam = sampler_func(width)
        end
    end
    v[:] = sam
    v
end

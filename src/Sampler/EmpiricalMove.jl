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

#validate(v::SamplerVariate{SliceTune{F}}) where {F<:SliceForm} =
#  validate(v, v.tune.width)
#
#validate(v::SamplerVariate{SliceTune{F}}, width::Int64) = v
#
#function validate(v::SamplerVariate{EmpiricalTune}, width::Int64)
#  n = length(v)
#  length(width) == n ||
#    throw(ArgumentError("length(width) differs from variate length $n"))
#  v
#end


#################### Sampler Constructor ####################

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
    v[:] = sam
    v
end

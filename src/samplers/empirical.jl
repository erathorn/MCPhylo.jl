#################### Empirical Move ####################

#################### Types and Constructors ####################


mutable struct EmpiricalTune{F<:Function,F2<:Function} <: SamplerTune
    logf::F
    samplerfun::F2
    width::Int64
    replacement::Bool

end

EmpiricalTune(
    x::Vector,
    logf::F,
    samplerfun::F2,
    width::Int64,
    replacement::Bool = false,
) where {F,F2} = EmpiricalTune{F,F2}(logf, samplerfun, width, replacement)

const EmpiricalVariate = Sampler{EmpiricalTune,T} where {T}


#################### Sampler Constructor ####################
"""

function Empirical(params::Symbol, width::Int64; replace::Bool=false)

This Sampler generates Samples from a distribution by selecting a random value
from this distribution. The move will always be accepted. Sampling can be done with or without replacement.
"""
function Empirical(params::Symbol, width::Int64; replace::Bool = false)

    tune = EmpiricalTune(Float64[], logpdf!, rand, width, replace)

    Sampler(Float64[], asvec(params), tune, Symbol[], false)
end


function sample!(
    v::Sampler{EmpiricalTune{F,F2},T},
    lpdf::Function;
    model::Model = model,
    kwargs...,
) where {T,F,F2}
    width = v.tune.width
    samfun = v.tune.samplerfun
    params = v.params
    sam = samfun(model, params, 1)
    if width > 1 && v.tune.replacement
        while length(Set(sam)) != width
            sam = samfun(model, params, 1)
        end
    end
    v[:] = sam
    v
end

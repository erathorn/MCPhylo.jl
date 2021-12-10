#################### Empirical Move ####################

#################### Types and Constructors ####################


mutable struct EmpiricalTune <: SamplerTune
    logf::Function
    samplerfun::Union{Function, Missing}
    width::Int64
    replacement::Bool

end

EmpiricalTune(x::Vector, logf::Function, samplerfun::Function, width::Int64, replacement::Bool=false)=
    EmpiricalTune(logf, samplerfun, width, replacement)

const EmpiricalVariate = Sampler{EmpiricalTune, T} where T


#################### Sampler Constructor ####################
"""

function Empirical(params::ElementOrVector{Symbol}, width::Int64; replace::Bool=false)

This Sampler generates Samples from a distribution by selecting a random value
from this distribution. The move will always be accepted. Sampling can be done with or without replacement.
"""
function Empirical(params::ElementOrVector{Symbol},
                width::Int64;
                replace::Bool=false)
    
    tune = EmpiricalTune(logpdf!, rand, width, replace)
    
    Sampler(Float64[], asvec(params), tune, Symbol[], false)
end


function sample!(v::EmpiricalVariate{T}, lpdf::Function; model::Model=model, kwargs...) where T
    width = v.tune.width
    samfun = v.tune.sampler_fun
    params = v.params
    sam = samfun(model, params, x)
    if width > 1 && v.tune.replacement
        while length(Set(sam)) != width
            sam = samfun(model, params, x)
        end
    end
    v[:] = sam
    v
end

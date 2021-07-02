
mutable struct PPHMCTune <: SamplerTune
    logfgrad::Union{Function,Missing}
    epsilon::Float64
    delta::Float64
    nleap::Int64
    moves::Vector{Int}

    PPHMCTune() = new()

    function PPHMCTune(x::Vector{T}, logfgrad::Union{Function,Missing}, epsilon::Float64, nleap::Int64, delta::Float64=0.003) where T <: GeneralNode
        new(logfgrad, epsilon, delta, nleap, Int[])
    end
end


#PPHMCTune(x::Vector{T}, logfgrad::Function, epsilon::Float64, nleap::Int64, delta::Float64=0.003) where T <: GeneralNode =
#    PPHMCTune(x, logfgrad, epsilon, nleap, delta)


const PPHMCVariate = SamplerVariate{PPHMCTune}



#################### Sampler Constructor ####################


"""
    PPHMC(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for PPHMC sampling. The Parameter is assumed to be
a tree.

Returns a `Sampler{PPHMCTune}` type object.

* params: stochastic node to be updated with the sampler.

* args...: additional keyword arguments to be passed to the PNUTSVariate constructor.
"""
function PPHMC(params::ElementOrVector{Symbol}, epsilon::Float64, nleap::Int64, delta::Float64)
    samplerfx = function (model::Model, block::Integer)
        block = SamplingBlock(model, block, true)

        f = let block = block
            (x, sz, ll, gr) -> mlogpdfgrad!(block, x, sz, ll, gr)
        end
        v = SamplerVariate(block, f, epsilon, nleap, delta)

        sample!(v::PPHMCVariate, f)

        relist(block, v)
    end
    Sampler(params, samplerfx, PPHMCTune())
end

sample!(v::PPHMCVariate; args...) = sample!(v, v.tune.logfgrad)

function sample!(v::PPHMCVariate, logfgrad::Function)
    mt = v.value[1]
    nl = size(mt)[1] - 1
    delta = v.tune.delta

    r = randn(nl)
    _, grad = logfgrad(mt, nl, true, true)
    ovnni = 0
    for i in 1:v.tune.nleap    
        mt, r, logf, grad, nni = refraction(mt, r, 1, grad, v.tune.epsilon, logfgrad, delta, nl)
        ovnni += nni
    end
    push!(v.tune.moves, ovnni)
    v.value[1] = mt
    v
end
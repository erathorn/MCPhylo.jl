
mutable struct PPHMCTune <: SamplerTune
    logfgrad::Union{Function,Missing}
    epsilon::Float64
    delta::Float64
    nleap::Int64
    m::Int64
    randomization::Bool
    adapter::Float64
    moves::Vector{Int}

    PPHMCTune() = new()

    function PPHMCTune(x::Vector{T}, logfgrad::Union{Function,Missing}, epsilon::Float64, nleap::Int64, delta::Float64=0.003) where T <: GeneralNode
        new(logfgrad, epsilon, delta, nleap, 0, true, 0.4,Int[])
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

        sample!(v::PPHMCVariate, f, model.iter <= model.burnin, model.burnin)

        relist(block, v)
    end
    Sampler(params, samplerfx, PPHMCTune())
end

sample!(v::PPHMCVariate; args...) = sample!(v, v.tune.logfgrad, bi)

function sample!(v::PPHMCVariate, logfgrad::Function, adapt::Bool, bi::Int64)
    mt = v.value[1]
    nl = size(mt)[1] - 1
    delta = v.tune.delta
    v.tune.m += 1
    r = randn(nl)
    blv = get_branchlength_vector(mt)
    set_branchlength_vector!(mt, molifier.(blv, delta))
    CU, grad = logfgrad(mt, nl, true, true)
    
    currH = CU - 0.5 * dot(r)
    ovnni = 0
    PU = 0.0
    fac = adapt ? (1-v.tune.adapter)^(1-v.tune.m*1/bi) : 1
    epsilon = v.tune.epsilon * fac
    nl = v.tune.randomization ? rand(1:v.tune.nleap) : v.tune.nleap
    for i in 1:nl
        mt, r, PU, grad, nni = refraction(mt, r, grad, epsilon, logfgrad, delta, nl)
        ovnni += nni
    end

    push!(v.tune.moves, ovnni)
    propH = PU - 0.5 * dot(r)
    
    if log(rand()) < propH -currH
        v.value[1] = mt
    end
    v
end
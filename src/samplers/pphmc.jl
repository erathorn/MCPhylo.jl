
mutable struct PPHMCTune{F<:Function,F2<:Function,G<:GradType} <: GradSampler{G}
    logf::F
    logfgrad::F2
    epsilon::Float64
    delta::Float64
    nleap::Int64
    m::Int64
    randomization::Bool
    adapter::Float64
    moves::Vector{Int}

    #PPHMCTune{F, F2}() where {F, F2} = new{F, F2}()

    function PPHMCTune(
        x::Vector{T},
        logf::F,
        logfgrad::F2,
        epsilon::Float64,
        nleap::Int64,
        ::Type{G};
        delta::Float64 = 0.003,
        adapter::Float64 = 0.4,
        randomization::Bool = true,
    ) where {T<:GeneralNode,F,F2,G<:GradType}
        new{F,F2,G}(logf, logfgrad, epsilon, delta, nleap, 0, randomization, adapter, Int[])
    end
end



const PPHMCVariate = Sampler{PPHMCTune,T} where {T<:GeneralNode}



#################### Sampler Constructor ####################


"""
    PPHMC(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for PPHMC sampling. The Parameter is assumed to be
a tree.

Returns a `Sampler{PPHMCTune}` type object.

* params: stochastic node to be updated with the sampler.

* args...: additional keyword arguments to be passed to the PNUTSVariate constructor.
"""
function PPHMC(
    params::ElementOrVector{Symbol},
    epsilon::Float64,
    nleap::Int64;
    delta::Float64 = 0.003,
    adapter::Float64 = 0.4,
    randomization::Bool = true,
)

    tune = PPHMCTune(
        GeneralNode[],
        logpdf!,
        logpdfgrad!,
        epsilon,
        nleap,
        provided,
        delta = delta,
        adapter = adapter,
        randomization = randomization,
    )
    Sampler(params, tune, Symbol[], false)
end


function sample!(
    v::Sampler{PPHMCTune{F,F2,G},Vector{T}},
    logfun::Function;
    grlpdf::Function,
    adapt::Bool,
    model::Model,
    kwargs...,
) where {F,F2,T<:GeneralNode,G}

    bi = model.burnin
    mt = v.value[1]
    nl = size(mt)[1] - 1
    delta = v.tune.delta
    v.tune.m += 1
    r = randn(nl)
    blv = get_branchlength_vector(mt)
    set_branchlength_vector!(mt, molifier.(blv, delta))
    lf, grad = grlpdf(mt)
    grad .*= scale_fac.(blv, delta)
    s = Tree_HMC_State(deepcopy(mt), r[:], grad[:], lf)
    currH = hamiltonian(s)
    ovnni = 0


    fac = adapt ? (1 - v.tune.adapter)^(1 - v.tune.m * 1 / bi) : 1
    epsilon = v.tune.epsilon * fac
    leaps = v.tune.randomization ? rand(1:v.tune.nleap) : v.tune.nleap
    att = 0
    for i = 1:leaps
        nni, att_nni = refraction!(s, epsilon, grlpdf, logfun, delta)
        ovnni += nni
        att += att_nni
    end

    push!(v.tune.moves, ovnni)

    propH = hamiltonian(s)
    if log(rand()) < propH - currH
        v.value[1] = deepcopy(s.x)
    end
    v
end
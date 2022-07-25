#################### No-U-Turn Sampler ####################

#################### Types and Constructors ####################

abstract type classic_nuts end
abstract type riemann_nuts end

const NUTS_Form = Union{classic_nuts, riemann_nuts}
mutable struct NUTSTune{N<:NUTS_Form, F<:Function} <: SamplerTune
  logf::F
  stepsizeadapter::NUTSstepadapter
  adapt::Bool
  epsilon::Float64
  tree_depth::Int
  acc_p_r::Vector{Int}

  NUTSTune{N, F}() where {N<:NUTS_Form, F<:Function} = new{N, F}()

  function NUTSTune{N}(x::Vector, epsilon::Real, logfgrad::F;
                    target::Real=0.75, tree_depth::Int=10) where {N<:NUTS_Form, F}
    new{N, F}(logfgrad,
    NUTSstepadapter(0,0,0,NUTS_StepParams(0.5,target,0.05,0.75,10,)),
    false,
    epsilon,
    tree_depth,
    Int[])
  end
end



#################### Sampler Constructor ####################

"""
    NUTS(params::ElementOrVector{Symbol}, classic_nuts; epsilon::Real = -Inf, kwargs...)

Construct a `Sampler` object for NUTS sampling, with the algorithm´s step size
parameter adaptively tuned during burn-in iterations. Parameters are assumed
to be continuous, but may be constrained or unconstrained. The `variant` parameter controls
whether the Riemannian variant of the NUTS algorithm (https://arxiv.org/abs/1304.1920) or the
classical variant (https://jmlr.org/papers/v15/hoffman14a.html) is used.

Returns a `Sampler{NUTSTune}` type object.

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* variant of the NUTS sampler. Options are
    * `classic_nuts` : classic version (default)
    * `riemann_nuts` : Riemannian version.

* `epsilon`: step size which can be presupplied. If the standard (-Inf) is chosen, it will be automatically adapted
"""
function NUTS(params::ElementOrVector{Symbol},::Type{F} = classic_nuts, epsilon::Real = -Inf, kwargs...) where F<:NUTS_Form
  tune = NUTSTune{F}(Float64[], epsilon, logpdfgrad!; kwargs...)
  Sampler(Float64[], params, tune, Symbol[], true)
end



#################### Sampling Functions ####################




"""
    sample!(v::NUTSVariate, logfgrad::Function; adapt::Bool=false)

Draw one sample from a target distribution using the NUTS sampler. Parameters
are assumed to be continuous and unconstrained.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::Sampler{NUTSTune{N, F}, T}, logfgrad::Function; adapt::Bool = false, kwargs...) where {T<: AbstractArray{<: Real}, N<:NUTS_Form, F<:Function}
    tune = v.tune
    adapter = tune.stepsizeadapter
    const_params = tune.stepsizeadapter.params
    if adapter.m == 0 && isinf(tune.epsilon)
        tune.epsilon = nutsepsilon(v.value, logfgrad, const_params.δ)
    end
    setadapt!(v, adapt)
    
    if tune.adapt
        adapter.m += 1
        
        nuts_sub!(v, tune.epsilon, logfgrad)

        adaptstat = adapter.metro_acc_prob > 1 ? 1 : adapter.metro_acc_prob
        
        adaptstat = const_params.δ - adaptstat
        
        η = 1.0/(adapter.m + const_params.t0)
        
        adapter.s_bar = (1.0 - η) * adapter.s_bar + η * adaptstat
        x = const_params.μ - adapter.s_bar * sqrt(adapter.m) / const_params.γ
        
        x_η = adapter.m^-const_params.κ
        adapter.x_bar = (1.0 - x_η) * adapter.x_bar + x_η * x
        tune.epsilon = exp(x)
        

    else
        if (adapter.m > 0)
            tune.epsilon = exp(adapter.x_bar)
        end

        nuts_sub!(v, tune.epsilon, logfgrad)
    end
    v
end


function setadapt!(v::Sampler{NUTSTune{N, F}, T}, adapt::Bool) where {T<: AbstractArray{<: Real}, N<:NUTS_Form, F<:Function}
    tune = v.tune
    if adapt && !tune.adapt
        tune.stepsizeadapter.m = 0
        tune.stepsizeadapter.params = update_step(tune.stepsizeadapter.params, log(10.0 * tune.epsilon))
    end
    tune.adapt = adapt
    v
end
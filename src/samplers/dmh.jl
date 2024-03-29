
mutable struct DMHTune{F<:Function, F2<:Function, F3<:Function} <: SamplerTune
    logf::F
    pseudolog::F2 # AuxLog
    condlike::F3
    datakeys::Vector{Symbol}
    m::Int64
    scale::Float64

    DMHTune{F, F2, F3}() where {F, F2, F3} = new{F, F2, F3}()

    function DMHTune(f::F, ps::F2, cl::F3, m::Int64, s::Float64) where {F, F2, F3}
        new{F, F2, F3}(f, ps, cl, Vector{Symbol}(), m, s)
    end
    function DMHTune(
        f::F,
        ps::F2,
        cl::F3,
        dk::Vector{Symbol},
        m::Int64,
        s::Float64,
    ) where {F, F2, F3}
        new{F, F2, F3}(f, ps, cl, dk, m, s)
    end
end

const DMHVariate = Sampler{DMHTune{F, F2, F3},T} where {F, F2, F3, T}

DMHTune(x::Vector, f, ps, cl, dk, m, s; args...) = DMHTune(f, ps, cl, dk, m, s; args...)

#################### Sampler Constructor ####################

function DMH(
    params::ElementOrVector{Symbol},
    m::Int64,
    scale::Float64 = 1.0,
    transform::Bool = true;
    args...,
)
    
    tune = DMHTune(logpdf!, pseudologpdf!, conditional_likelihood!, m, scale)

    Sampler(params, tune, Symbol[], transform)
end



function sample!(v::Sampler{DMHTune{F, F2, F3},T}, lpdf::Function; model::Model, kwargs...) where {T, F, F2, F3}
    # 1. propose theta prime
    # 2. Generate the auxiliary variable using theta prime

    tune = v.tune
    
    targets = keys(model, :target, v.params)
    data = unlist(model[targets[1]])
    data = reshape(data, size(getindex(model, targets[1])))
    nfeatures, nlang = size(data)
    
    @assert v.tune.m >= nlang

    # store old value for possible future reference
    θ = deepcopy(v.value)
    # this way of generating theta_prime from the current values of theta
    # takes care of the transition probability from theta_prime to theta and vice versa
    # the values equal and will cancel out.
    θ_prime = θ + v.tune.scale .* rand(Normal(0.0, 1.0), length(v))

    
    # calculate logpdf values
    lf_xt = lpdf(v.value)

    
    lf_xtp =lpdf(θ_prime)

    # prematurely assign, to allow computations to go through properly
    v[:] = θ_prime

    # Sample new pseudo observations based on the θ_prime
    y = inner_sampler(v, deepcopy(data), v.params, targets, model)

    # calculate logpdfs based on pseudo observations
    lf_ytp = tune.pseudolog(model, v.value, y, v.params, targets, true)
    lf_yt = tune.pseudolog(model, θ, y, v.params, targets, true)

    # calculate acceptance probability (proposal distribution?)
    r = exp((lf_yt + lf_xtp) - (lf_xt + lf_ytp))

    # RWM acceptance logic is a bit reversed here. Above the new value is
    # prematurely assigned to allow computations in the inner sampler to go through.
    # If the sample is accepted nothing needs to be done anymore, otherwise the
    # old value will be reassigned.

    if rand() > r
        # sample is rejected, so use the old value.
        v[:] = θ
    end
    v
end

function inner_sampler(
    v::Sampler{DMHTune{F, F2, F3},T},
    X::Array{N,2},
    params,
    targets,
    model,
)::Array{N,2} where {N<:Real, F, F2, F3, T}
    nfeatures, nlangs = size(X)
    counter = zero(Int64)

    while true
        random_language_idx = shuffle(1:nlangs)
        random_feature_idx = shuffle(1:nfeatures)
        @inbounds for i in random_language_idx, f in random_feature_idx

            probs = v.tune.condlike(model, X, params, targets, i, f)
            new_x = sample(1:length(probs), Weights(probs))
            X[f, i] = new_x
        
            counter += 1
            if counter > v.tune.m
                return X
            end
        end
    end
end

#################### Missing Values Sampler ####################

#################### Types and Constructors ####################

struct MISSTune{F<:Function} <: SamplerTune
    logf::F
    dims::Tuple
    valueinds::AbstractArray
    distrinds::AbstractArray

end

function MISSTune(f::F) where {F}
    MISSTune(f, (), Int[], Int[])
end


function MISSTune(s::AbstractStochastic, gm::Function)
    MISSTune(s.distr, s.value, gm)
end

function MISSTune(d::Distribution, v, gm::Function)
    MISSTune(gm, (), findall(isnan.(v)), Int[])
end

function MISSTune(D::Array{UnivariateDistribution}, v::Array, gm::Function)
    inds = findall(isnan.(v))
    MISSTune(gm, dims(D), inds, inds)
end

function MISSTune(D::Array{MultivariateDistribution}, v::Array, gm::Function)
    isvalueinds = falses(size(v))
    isdistrinds = falses(size(D))
    for sub in CartesianIndices(size(D))
        n = length(D[sub])
        for i = 1:n
            if isnan(v[sub, i])
                isvalueinds[sub, i] = isdistrinds[sub] = true
            end
        end
    end
    MISSTune(gm, dims(D), findall(isvalueinds), findall(isdistrinds))
end

const MISSVariate = Sampler{MISSTune,T} where {T}

#################### Sampler Constructor ####################
"""
    MISS(params::ElementOrVector{Symbol})

Construct a `Sampler` object to sampling missing output values. The constructor
should only be used to sample stochastic nodes upon which no other stochastic
node depends. So-called ‘output nodes’ can be identified with the `keys()`
function. Moreover, when the `MISS` constructor is included in a vector of
`Sampler` objects to define a sampling scheme, it should be positioned at the
beginning of the vector. This ensures that missing output values are updated
before any other samplers are executed.

Returns a `Sampler{Dict{Symbol, MISSTune}}` type object.

* `params`: stochastic node(s) that contain missing values (`NaN`) to be updated with the sampler.
"""
function MISS(params::Symbol)
    lf(m, args...) = m
    params = asvec(params)

    Sampler(Float64[], params, MISSTune(lf), Symbol[], false)
end


#################### Sampling Functions ####################


function sample!(
    v::Sampler{MISSTune{F},T},
    get_model::Function;
    gen::Int,
    model::Model,
    kwargs...,
) where {T,F}
    params = v.params
    if gen == 1
        for key in params
            v.tune = MISSTune(model[key], get_model)
        end
    end
    for key in params
        node = model[key]

        node[v.tune.valueinds] = rand(node, v.tune)
        update!(model, node.targets)
    end
    nothing
end



rand(s::AbstractStochastic, miss::MISSTune) = rand_sub(s.distr, miss)

function rand_sub(d::Distribution, miss::MISSTune)
    x = rand(d)
    Float64[x[i] for i in miss.valueinds]
end

function rand_sub(D::Array{UnivariateDistribution}, miss::MISSTune)
    Float64[rand(d) for d in D[miss.distrinds]]
end

function rand_sub(D::Array{MultivariateDistribution}, miss::MISSTune)
    X = Array{Float64}(undef, miss.dims)
    for i in miss.distrinds
        d = D[i]
        X[ind2sub(D, i)..., 1:length(d)] = rand(d)
    end
    X[miss.valueinds]
end

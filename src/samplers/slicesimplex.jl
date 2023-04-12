#################### Slice Simplex Sampler ####################

#################### Types and Constructors ####################

mutable struct SliceSimplexTune{F<:Function} <: SamplerTune
    logf::F
    scale::Float64

    SliceSimplexTune(x::Vector, logf::F; scale::Real = 1.0) where {F} = new{F}(logf, scale)
end

SliceSimplexTune(x::Vector; args...) = SliceSimplexTune(x, identity; args...)

const SliceSimplexVariate = Sampler{SliceSimplexTune,T} where {T}

function validate(v::Sampler{SliceSimplexTune{F},T}) where {F,T}
    0 < v.tune.scale <= 1 || throw(ArgumentError("scale is not in (0, 1]"))
    validatesimplex(v)
end


#################### Sampler Constructor ####################
"""
    SliceSimplex(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for which SliceSimplex sampling is to be applied
separately to each of the supplied parameters. Parameters are assumed to be
continuous and constrained to a simplex.

Returns a `Sampler{SliceSimplexTune}` type object.

* `params`: stochastic node(s) to be updated with the sampler.

* `args...`: additional keyword arguments to be passed to the `SliceSimplexVariate` constructor.
"""
function SliceSimplex(params::Symbol; scale::Real = 1.0, args...)
    params = asvec(params)
    tune = SliceSimplexTune(Float64[], logpdf!, scale = scale)
    Sampler(params, tune, Symbol[], false)
end




function SliceSimplex_sub!(d::MultivariateDistribution, sim::Function, logf::Function)
    inds = 1:length(d)
    sim(inds, v -> logf(d, v, inds))
end

function SliceSimplex_sub!(
    D::Array{MultivariateDistribution},
    sim::Function,
    logf::Function,
)
    inds = 0:0
    @inbounds for i = 1:length(D)
        d = D[i]
        inds = last(inds) .+ (1:length(d))
        sim(inds, v -> logf(d, v, inds))
    end
end

function SliceSimplex_sub!(d, sim::Function, logf::Function)
    throw(ArgumentError("unsupported distribution structure $(typeof(d))"))
end


#################### Sampling Functions ####################

#sample!(v::SliceSimplexVariate) = sample!(v, v.tune.logf)
"""
    sample!(v::SliceSimplexVariate, logf::Function)

Draw one sample from a target distribution using the SliceSimplex sampler.
Parameters are assumed to be continuous and constrained to a simplex.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::Sampler{SliceSimplexTune{F},T}, logf::Function; args...) where {F,T}


    p0 = logf(v.value) + log(rand())

    d = Dirichlet(length(v.value), 1)
    vertices = makefirstsimplex(v.value, v.tune.scale)
    vb = vertices \ v.value
    xb = rand(d)
    x = vertices * xb
    while any(x .< 0.0) || any(x .> 1.0) || logf(x) < p0
        shrinksimplex!(vertices, vb, xb, v.value, x)
        vb = vertices \ v.value
        xb = rand(d)
        x = vertices * xb
    end
    v[:] = x

    v
end


function makefirstsimplex(x::AbstractVector{T}, scale::Real)::Matrix{T} where {T}
    vertices = Matrix{T}(I, length(x), length(x))
    vertices[:, 2:end] += (1.0 - scale) * (vertices[:, 1] .- vertices[:, 2:end])
    vertices .+ x .- vertices * rand(Dirichlet(fill!(similar(x), 1)))
end


function shrinksimplex!(
    vertices::S,
    bx::T,
    bc::T,
    cx::T,
    cc::T,
) where {T<:AbstractVector{<:Real},S<:AbstractArray{<:Real}}
    @inbounds @fastmath for i in findall(bc .< bx)
        inds = [1:(i-1); (i+1):size(vertices, 2)]
        vertices[:, inds] += bc[i] * (vertices[:, i] .- vertices[:, inds])
        bc = vertices \ cc
    end

    vertices
end

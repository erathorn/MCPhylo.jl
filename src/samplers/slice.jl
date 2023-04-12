#################### Slice Sampler ####################

#################### Types and Constructors ####################

const SliceForm = Union{Univariate,Multivariate}

mutable struct SliceTune{F<:SliceForm} <: SamplerTune
    logf::Union{Function,Missing}
    width::Union{Float64,Vector{Float64}}

    SliceTune{F}() where {F<:SliceForm} = new{F}()

    SliceTune{F}(x::Vector, width) where {F<:SliceForm} = SliceTune{F}(x, width, missing)

    SliceTune{F}(
        x::Vector,
        width::Real,
        logf::Union{Function,Missing},
    ) where {F<:SliceForm} = new{F}(logf, Float64(width))

    SliceTune{F}(
        x::Vector,
        width::Vector,
        logf::Union{Function,Missing},
    ) where {F<:SliceForm} = new{F}(logf, convert(Vector{Float64}, width))
    SliceTune{F}(width, logf::Union{Function,Missing}) where {F<:SliceForm} =
        new{F}(logf, width)
end


const SliceUnivariate = Sampler{SliceTune{Univariate},R} where {R}
const SliceMultivariate = Sampler{SliceTune{Multivariate},R} where {R}

validate(v::Sampler{SliceTune{F}}) where {F<:SliceForm} = validate(v, v.tune.width)

validate(v::Sampler{SliceTune{F}}, width::Float64) where {F<:SliceForm} = v

function validate(v::Sampler{SliceTune{F}}, width::Vector) where {F<:SliceForm}
    n = length(v)
    length(width) == n ||
        throw(ArgumentError("length(width) differs from variate length $n"))
    v
end


#################### Sampler Constructor ####################
"""
    Slice(params::ElementOrVector{Symbol},
                width::ElementOrVector{T},
                ::Type{F}=Multivariate;
                transform::Bool=false) where {T<:Real, F<:SliceForm}

Construct a `Sampler` object for Slice sampling. Parameters are assumed to be
continuous, but may be constrained or unconstrained.

Returns a `Sampler{SliceTune{Univariate}}` or `Sampler{SliceTune{Multivariate}}`
type object if sampling univariately or multivariately, respectively.

  * `params`: stochastic node(s) to be updated with the sampler.

  * `width`: scaling value or vector of the same length as the combined elements of nodes `params`, defining initial widths of a hyperrectangle from which to simulate values.

  * `F` : sampler type. Options are
      * `:Univariate` : sequential univariate sampling of parameters.
      * `:Multivariate` : joint multivariate sampling.

  * `transform`: whether to sample parameters on the link-transformed scale (unconstrained parameter space). If `true`, then constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function, and `width` is interpreted as being relative to the unconstrained parameter space. Otherwise, sampling is relative to the untransformed space.
"""
function Slice(
    params::ElementOrVector{Symbol},
    width::ElementOrVector{T},
    ::Type{F} = Multivariate;
    transform::Bool = false,
) where {T<:Real,F<:SliceForm}
    tune = SliceTune{F}(width, logpdf!)
    Sampler(params, tune, Symbol[], transform)
end


#################### Sampling Functions ####################

#sample!(v::Union{SliceUnivariate,SliceMultivariate}, logf::Function; args...) = sample!(v, logf)
"""
    sample!(v::Union{SliceUnivariate, SliceMultivariate}, logf::Function)

Draw one sample from a target distribution using the Slice univariate or
multivariate sampler. Parameters are assumed to be continuous, but may be
constrained or unconstrained.

Returns `v` updated with simulated values and associated tuning parameters.
"""

function sample!(
    v::SliceUnivariate{Vector{T}},
    logf::Function;
    kwargs...,
) where {T<:GeneralNode}
    tree = v.value[1]

    logf0 = logf(v.value)
    blv = get_branchlength_vector(tree)

    n = length(blv)
    lower = blv - v.tune.width .* rand(n)
    lower[lower.<0.0] .= 0.0
    upper = lower .+ v.tune.width

    @inbounds for i = 1:n
        p0 = logf0 + log(rand())

        x = blv[i]
        blv[i] = rand(Uniform(lower[i], upper[i]))
        while true
            set_branchlength_vector!(tree, blv)
            logf0 = logf([tree])
            logf0 < p0 || break
            value = blv[i]
            if value < x
                lower[i] = value
            else
                upper[i] = value
            end
            blv[i] = rand(Uniform(lower[i], upper[i]))
        end
    end

    v
end


function sample!(
    v::SliceMultivariate{Vector{T}},
    logf::Function;
    kwargs...,
) where {T<:GeneralNode}
    tree = v.value[1]

    p0 = logf(v.value) + log(rand())
    blv = get_branchlength_vector(tree)
    org = deepcopy(blv)


    n = length(blv)
    lower = blv - v.tune.width .* rand(n)
    lower[lower.<0.0] .= 0.0
    upper = lower .+ v.tune.width

    blv = v.tune.width .* rand(n) + lower
    set_branchlength_vector!(tree, blv)
    @inbounds while logf([tree]) < p0
        for i = 1:n
            value = blv[i]
            if value < org[i]
                lower[i] = value
            else
                upper[i] = value
            end
            blv[i] = rand(Uniform(lower[i], upper[i]))
        end
        set_branchlength_vector!(tree, blv)
    end

    v.value[1] = tree
    v
end



function sample!(v::SliceUnivariate{Vector{Float64}}, logf::Function; kwargs...)
    logf0 = logf(v.value)

    n = length(v)
    lower = v - v.tune.width .* rand(n)
    upper = lower .+ v.tune.width

    @inbounds for i = 1:n
        p0 = logf0 + log(rand())

        x = v[i]
        v[i] = rand(Uniform(lower[i], upper[i]))
        while true
            logf0 = logf(v.value)
            logf0 < p0 || break
            value = v[i]
            if value < x
                lower[i] = value
            else
                upper[i] = value
            end
            v[i] = rand(Uniform(lower[i], upper[i]))
        end
    end

    v
end


function sample!(v::SliceMultivariate{Vector{Float64}}, logf::Function; kwargs...)
    p0 = logf(v.value) + log(rand())

    n = length(v)
    lower = v - v.tune.width .* rand(n)
    upper = lower .+ v.tune.width

    x = v.tune.width .* rand(n) + lower
    while logf(x) < p0
        @inbounds for i = 1:n
            value = x[i]
            if value < v.value[i]
                lower[i] = value
            else
                upper[i] = value
            end
            x[i] = rand(Uniform(lower[i], upper[i]))
        end
    end
    v[:] = x

    v
end

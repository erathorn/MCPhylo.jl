#################### Sampler ####################

const samplerfxargs = [(:model, MCPhylo.Model), (:block, Integer)]

const chunksize = 40

#################### Types and Constructors ####################

struct NullFunction end

struct SamplingBlock
    model::Model
    index::Int
    transform::Bool

    SamplingBlock(model::Model, index::Integer = 0, transform::Bool = false) =
        new(model, index, transform)
end


Sampler(param::Symbol, args...) = Sampler([param], args...)
"""
    Sampler(params::Vector{Symbol}, f::Function, tune::Any=Dict())
Construct a `Sampler` object that defines a sampling function for a block of
stochastic nodes.

Returns a `Sampler{typeof(tune)}` type object.

* `params`: node(s) being block-updated by the sampler.

* `f`: function for the `eval` field of the constructed sampler and whose arguments are the other model nodes upon which the sampler depends, typed argument `model::Model` that contains all model nodes, and/or typed argument `block::Integer` that is an index identifying the corresponding sampling function in a vector of all samplers for the associated model. Through the arguments, all model nodes and fields can be accessed in the body of the function. The function may return an updated sample for the nodes identified in its `params` field. Such a return value can be a structure of the same type as the node if the block consists of only one node, or a dictionary of node structures with keys equal to the block node symbols if one or more. Alternatively, a value of `nothing` may be returned. Return values that are not `nothing` will be used to automatically update the node values and propagate them to dependent nodes. No automatic updating will be done if `nothing` is returned.

* `tune`: tuning parameters needed by the sampling function.
"""
function Sampler(params::Vector{Symbol}, f::Function, tune::Any = Dict())
    Sampler(params, modelfx(samplerfxargs, f), tune, Symbol[])
end

function Sampler(
    params::Vector{Symbol},
    tune::SamplerTune,
    targets::Vector{Symbol},
    transform::Bool = false,
)
    Sampler(Float64[], params, tune, targets, transform)
end

function Sampler(v::T, s::Sampler{R,X})::Sampler{R,T} where {R<:SamplerTune,X,T}
    Sampler(v, s.params, s.tune, s.targets, s.transform)
end


#################### Base Methods ####################

function Base.show(io::IO, s::Sampler)
    print(io, "An object of type \"$(summary(s))\"\n")
    print(io, "Sampling Block Nodes:\n")
    show(io, s.params)
    #print(io, "\n\n")
    #show(io, "text/plain", first(code_typed(s.eval)))
    println(io)
end

function showall(io::IO, s::Sampler)
    show(io, s)
    print(io, "\nTuning Parameters:\n")
    show(io, s.tune)
    print(io, "\n\nTarget Nodes:\n")
    show(io, s.targets)
end


#################### Variate Validators ####################

validate(v::Sampler) = v

function validatebinary(v::Sampler)
    all(insupport(Bernoulli, v)) || throw(ArgumentError("variate is not a binary vector"))
    v
end

function validatesimplex(v::Sampler)
    isprobvec(v) || throw(ArgumentError("variate is not a probability vector $v"))
    v
end


#################### Simulation Methods ####################

function gradlogpdf!(
    block::SamplingBlock,
    x::AbstractArray{T},
    dtype::Symbol = :forward,
) where {T<:Real}
    gradlogpdf!(block.model, x, block.index, block.transform, dtype = dtype)
end


function gradlogpdf!(block::SamplingBlock, x::N) where {N<:GeneralNode}
    gradlogpdf!(block.model, x, block.index, block.transform)
end

function logpdf!(block::SamplingBlock, x::AbstractArray{T}) where {T<:Real}
    logpdf!(block.model, x, block.index, block.transform)
end

function pseudologpdf!(
    block::SamplingBlock,
    x::AbstractArray{T},
    y::AbstractArray,
) where {T<:Real}
    pseudologpdf!(block.model, x, y, block.index, block.transform)
end

function conditional_likelihood!(
    block::SamplingBlock,
    x::AbstractArray{T},
    args...,
) where {T<:Real}
    conditional_likelihood!(block.model, x, block.index, args...)
end

function logpdf!(block::SamplingBlock, x::AbstractArray{T}) where {T<:GeneralNode}
    logpdf!(block.model, x[1])
end

function logpdf!(block::SamplingBlock, x::T) where {T<:GeneralNode}
    logpdf!(block.model, x, block.index, block.transform)
end

function rand!(block::SamplingBlock, x::Int64)
    rand!(block.model, x, block.index)
end



function _gradlogpdf!(m::Model, x::AbstractArray, block::Integer, dtype::Symbol = :provided)
    targets = keys(m, :target, block)

    node = m[targets[1]]
    update!(node, m)
    return gradlogpdf(node)

end

function logpdfgrad!(
    block::SamplingBlock,
    x::AbstractVector{T},
    dtype::Symbol,
) where {T<:Real}
    grad = gradlogpdf!(block, x, dtype)
    logf = logpdf!(block, x)
    (logf, ifelse.(isfinite.(grad), grad, 0.0))
end


function logpdfgrad!(
    m::Model,
    x::AbstractVector{T},
    params::ElementOrVector{Symbol},
    target::ElementOrVector{Symbol},
    transform::Bool,
) where {T<:Real}
    ll::Float64 = 0.0

    function lf(y)
        let m1 = m, para = params, tar = target, tr = transform
            lp = logpdf!(m1, y, para, tar, tr)
            ll = ForwardDiff.value(lp)
            lp
        end
    end
    #ll1 = lf(x)
    #grad1 = FiniteDiff.finite_difference_gradient(lf, x)
    chunk = ForwardDiff.Chunk(min(length(x), chunksize))
    config = ForwardDiff.GradientConfig(lf, x, chunk)
    grad = ForwardDiff.gradient!(similar(x), lf, x, config)

    #(ll, ifelse.(isfinite.(grad), grad, 0.0))
    ll, grad
end



function logpdfgrad!(
    m::Model,
    x::T,
    params::ElementOrVector{Symbol},
    target::ElementOrVector{Symbol},
    transform::Bool,
) where {T<:GeneralNode}

    m[params] = relist(m, x, params, transform)

    # likelihood
    v, grad = gradlogpdf(m, target)
    # prior
    vp, gradp = gradlogpdf(m[params[1]], x)

    vp + v, gradp .+ grad
end


function logpdf!(
    m::Model,
    x::T,
    params::ElementOrVector{Symbol},
    target::ElementOrVector{Symbol},
    transform::Bool,
) where {T<:GeneralNode}

    m[params] = relist(m, x, params, transform)

    # likelihood
    v = logpdf(m, target)
    # prior
    vp = logpdf(m[params[1]], x)

    vp + v
end





#################### unlist and relist functionality ####################

function unlist(block::SamplingBlock)
    unlist(block.model, block.index, block.transform)
end

function relist(block::SamplingBlock, x::AbstractVariate) where {T<:Real}
    relist(block.model, x.value, block.index, block.transform)
end


#################### Auxiliary Functions ####################

asvec(x::Union{Number,Symbol}) = [x]
asvec(x::Vector) = x

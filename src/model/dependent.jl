#################### Dependent ####################

const depfxargs = [(:model, MCPhylo.Model)]


#################### Base Methods ####################


dims(d::AbstractDependent) = size(d)

function names(d::AbstractDependent)
    names(d, d.symbol)
end
"""
    setmonitor!(d::AbstractDependent, monitor::Bool)

Specify node elements to be included in monitored MCMC sampler output.

Returns d with its monitor field updated to reflect the specified monitoring.

* `d` : node whose elements contain sampled MCMC values.

* `monitor` : boolean indicating whether all elements are monitored.

"""
function setmonitor!(d::AbstractDependent, monitor::Bool)
    value = monitor ? Int[0] : Int[]
    setmonitor!(d, value)
end
"""
    setmonitor!(d::AbstractDependent, monitor::Vector{Int})

Specify node elements to be included in monitored MCMC sampler output.

Returns d with its monitor field updated to reflect the specified monitoring.

* `d` : node whose elements contain sampled MCMC values.

* `monitor` : vector of element-wise indices of elements to monitor.
"""
function setmonitor!(d::AbstractDependent, monitor::Vector{Int})
    values = monitor
    if !isempty(monitor)
        n = isa(d, TreeVariate) ? length(unlist_tree(d)) : length(unlist(d))
        if n > 0
            if monitor[1] == 0
                values = collect(1:n)
            elseif minimum(monitor) < 1 || maximum(monitor) > n
                throw(BoundsError())
            end
        end
    end
    d.monitor = values
    d
end
#
"""
    setmonitor!(d::TreeVariate, monitor::Bool)

Specify node elements to be included in monitored MCMC sampler output.

Returns `d` with its `monitor` field updated to reflect the specified monitoring.

* `d` : node whose elements contain sampled MCMC values.

* `monitor` : boolean indicating whether all elements are monitored.
"""
function setmonitor!(d::TreeVariate, monitor::Bool)
    value = monitor ? Int[0] : Int[]
    setmonitor!(d, value)
end

"""
    setmonitor!(d::TreeVariate, monitor::Vector{Int})

Specify node elements to be included in monitored MCMC sampler output.

Returns `d` with its `monitor` field updated to reflect the specified monitoring.

* `d` : node whose elements contain sampled MCMC values.

* `monitor` : vector of element-wise indices of elements to monitor.
"""
function setmonitor!(d::TreeVariate, monitor::Vector{Int})
    values = monitor
    if !isempty(monitor)
        n = length(unlist_tree(d))
        if n > 0
            if monitor[1] == 0
                values = collect(1:n)
            elseif minimum(monitor) < 1 || maximum(monitor) > n
                throw(BoundsError())
            end
        end
    end
    d.monitor = values
    d
end



#################### Distribution Fallbacks ####################

unlist(d::AbstractDependent, transform::Bool = false) = unlist(d, d.value, transform)

unlist(d::AbstractDependent, x::Real, transform::Bool = false) = [x]

unlist(d::AbstractDependent, x::AbstractArray, transform::Bool = false) = vec(x)

relist(d::AbstractDependent, x::AbstractArray, transform::Bool = false) =
    relistlength(d, x, transform)[1]

logpdf(d::AbstractDependent, x=nothing , transform::Bool=false) = 0.0

gradlogpdf(d::AbstractDependent, x=nothing, transform::Bool=false) = 0.0



#################### Logical ####################

@promote_scalarvariate Logical{<:Real}


#################### Constructors ####################

"""
    Logical(f::Function, monitor::Union{Bool, Vector{Int}}=true)

Constructor for a Logical model node. This function assumes the output of the
logical operation to be a scalar.

* `f`: Function specifying the deterministic operation performed on its arguments. These arguments are other nodes of the model.

* `monitor` : Indicates whether the results should be monitored, i.e.saved.
"""
function Logical(f::Function, monitor::Union{Bool,Vector{Int}} = true)
    value = Float64(NaN)
    fx, src = modelfxsrc(depfxargs, f)
    l = Logical(value, :nothing, Int[], fx, src, Symbol[])
    setmonitor!(l, monitor)
end

Logical(f::Function, d::Integer, args...) = Logical(d, f, args...)

"""
    Logical(d::Integer, f::Function, monitor::Union{Bool, Vector{Int}}=true)

Constructor for a Logical model node.

* `d` : Specifies the dimension of the output.

* `f` : Specifies the deterministic operation performed on its arguments. These arguments are other nodes of the model.

* `monitor` : Indicates whether the results should be monitored, i.e.saved.
"""
function Logical(d::Integer, f::Function, monitor::Union{Bool,Vector{Int}} = true)
    value = Array{Float64}(undef, fill(0, d)...)
    fx, src = modelfxsrc(depfxargs, f)
    l = Logical(value, :nothing, Int[], fx, src, Symbol[])
    setmonitor!(l, monitor)
end

Logical(f::Function, d::T, args...) where {T<:GeneralNode} = Logical(d, f, args...)

"""
    Logical(d::T, f::Function, monitor::Union{Bool, Vector{Int}}=true) where T<:GeneralNode

Constructor for a Logical model node, which can hold a Node structure, i.e. a tree.

* `f` is a function specifying the deterministic operation performed on its arguments. These arguments are other nodes of the model.

* `monitor` indicates whether the results should be monitored, i.e. saved.
"""
function Logical(
    d::T,
    f::Function,
    monitor::Union{Bool,Vector{Int}} = true,
) where {T<:GeneralNode}
    value = Node()
    fx, src = modelfxsrc(depfxargs, f)
    l = Logical(value, :nothing, Int[], fx, src, Symbol[])
    setmonitor!(l, monitor)
end

ScalarLogical(x::T) where {T<:Real} = x

function Logical(a::Logical, value::T)::Logical{T} where T
    Logical(value, a.symbol, a.monitor, a.eval, a.sources, a.targets)
end

#################### Updating ####################
"""
    setinits!(l::AbstractLogical, m::Model, ::Any=nothing)

Set initial values for a logical node.

Returns the result of a call to setmonitor!(l, l.monitor) or setmonitor!(d, d.monitor).

* `l` : logical node to which to assign initial values.

* `m` : model containing the node.
"""
function setinits(l::Logical, m::Model, ::Any = nothing)
    l.value = l.eval(m)
    setmonitor!(l, l.monitor)
    l
end

"""
    update!(l::AbstractLogical, m::Model)

Update the values of a logical node according to its relationship with others in a model.

Returns the node with its values updated.

* `l` : logical node to update.

* `m` : model containing the node.
"""
function update!(l::T, m::Model) where T <: AbstractLogical   
    l1 = Logical(l.eval(m), l.symbol, l.monitor, l.eval, l.sources, l.targets)
    l1
end


#################### Distribution Methods ####################

relistlength(d::Logical{<:Real}, x::AbstractArray, transform::Bool = false) = (x[1], 1)

function relistlength(d::Logical{<:A}, x::A, transform::Bool = false) where  A<:AbstractArray
    n = length(d)
    value = reshape(x[1:n], size(d))
    (value, n)
end


#################### Stochastic ####################

#################### Base Methods ####################

@promote_scalarvariate Stochastic{<:Real}

function showall(io::IO, s::AbstractStochastic)
    show(io, s)
    print(io, "\n\nDistribution:\n")
    show(io, s.distr)
    print(io, "\nFunction:\n")
    show(io, "text/plain", first(code_typed(s.eval)))
    print(io, "\n\nSource Nodes:\n")
    show(io, s.sources)
    print(io, "\n\nTarget Nodes:\n")
    show(io, s.targets)
end


#################### Constructors ####################

"""
    Stochastic(f::Function, monitor::Union{Bool, Vector{Int}}=true)

Constructor for a Stochastic model node. This function assumes the output of the
logical operation to be scalar.

* `f` : Specifies the distributional relationship between the arguments and the node. These arguments are other nodes of the model.

* `monitor` : Indicates whether the results should be monitored, i.e. saved.
"""
function Stochastic(f::Function, monitor::Union{Bool,Vector{Int}} = true)
    value = Float64(NaN)
    fx, src = modelfxsrc(depfxargs, f)
    s = Stochastic(
        value,
        :nothing,
        Int[],
        fx,
        src,
        Symbol[],
        NullUnivariateDistribution(),
        -Inf
    )
    setmonitor!(s, monitor)
end

Stochastic(f::Function, d::Integer, args...) = Stochastic(d, f, args...)


function Stochastic(d::Integer, f::Function, monitor::Union{Bool,Vector{Int}} = true)

    value = Array{Float64}(undef, fill(0, d)...)

    fx, src = modelfxsrc(depfxargs, f)
    s = Stochastic(
        value,
        :nothing,
        Int[],
        fx,
        src,
        Symbol[],
        NullUnivariateDistribution(),
        -Inf
    )
    setmonitor!(s, monitor)
end


Stochastic(f::Function, d::T, args...) where {T<:GeneralNode} = Stochastic(d, f, args...)


"""

    Stochastic(d::N, f::Function, monitor::Union{Bool, Vector{Int}}=true) where N<:GeneralNode

Constructor for a Stochastic model node, which can hold a Node structure, i.e. a tree.

* `f` : Specifies the distributional reslationship between the arguments and the node.
These arguments are other nodes of the model.

* `monitor` : Indicates whether the results should be monitored, i.e. saved.
"""
function Stochastic(
    d::N,
    f::Function,
    monitor::Union{Bool,Vector{Int}} = true,
) where {N<:GeneralNode}
    value = Node()
    fx, src = modelfxsrc(depfxargs, f)
    s = Stochastic(
        value,
        :nothing,
        Int[],
        fx,
        src,
        Symbol[],
        NullUnivariateDistribution(),
        -Inf,
    )
    setmonitor!(s, monitor)
end

function Stochastic(a::Stochastic{S, F, D}, value::T)::Stochastic{T, F, D} where {T<:Union{Real,AbstractArray{T1,N} where {T1<:Real,N},GeneralNode}, F<:Function, S<:Union{Real,AbstractArray{T1,N} where {T1<:Real,N},GeneralNode}, D<:DistributionStruct}
    Stochastic(value, a.symbol, a.monitor, a.eval, a.sources, a.targets, a.distr, a.lpdf)
end

function Stochastic(a::Stochastic{T, F, S}, distr::D)::Stochastic{T, F, D} where {T<:Union{Real,AbstractArray{T1,N} where {T1<:Real,N},GeneralNode}, F<:Function, S<:DistributionStruct, D<:DistributionStruct}
    Stochastic(a.value, a.symbol, a.monitor, a.eval, a.sources, a.targets, distr, a.lpdf)
end
#################### Updating ####################
"""
    setinits!(s::ScalarStochastic, m::Model, x::Real)

Set initial values for a stochastic node.

Returns the node with its assigned initial values.

* `s` : stochastic node to which to assign initial values.

* `m` : model containing the node.

* `x` : values to assign to the node.
"""
function setinits(s::Stochastic{T}, m::Model, x::R) where R <: Real where T
    s.value = convert(Float64, x)
    s = Stochastic(s, s.eval(m))
    setmonitor!(s, s.monitor)
    s
end
"""
    setinits!(s::ArrayStochastic, m::Model, x::DenseArray)

Set initial values for a stochastic node.

Returns the node with its assigned initial values.

* `s` : stochastic node to which to assign initial values.

* `m` : model containing the node.

* `x` : values to assign to the node.
"""
function setinits(s::Stochastic{<:DenseArray}, m::Model, x::DenseArray)
  s.value = convert(typeof(s.value), copy(x))
  s = Stochastic(s, s.eval(m))
  if isa(s.distr, PhylogeneticDistribution)
    distrdims = dims(s.distr)
    for (ind, di) in enumerate(dims(s))
      if ind != 2
        di != distrdims[ind] && throw(DimensionMismatch("incompatible distribution for stochastic node"))
      end
    end
  elseif !isa(s.distr, UnivariateDistribution) && dims(s) != dims(s.distr)
    throw(DimensionMismatch("incompatible distribution for stochastic node $(s.symbol). Expected $(dims(s.distr)), got$(dims(s))."))
  end
  setmonitor!(s, s.monitor)
  return s
end


"""
    setinits!(d::TreeStochastic, m::Model, x::T) where {T<:GeneralNode}

Set initial values for a stochastic node.

Returns the node with its assigned initial values.

* `s` or `d` : stochastic node to which to assign initial values.

* `m` : model containing the node.

* `x` : values to assign to the node.
"""
function setinits(d::Stochastic{T}, m::Model, x::T) where {T<:GeneralNode}
    d.value = deepcopy(x)
    d = Stochastic(d, d.eval(m))
    insupport(d.distr, x) || throw(
        ArgumentError("The supplied tree does not match the topological tree constraints."),
    )
    setmonitor!(d, d.monitor)
    d
end # function


"""
    update!(s::AbstractStochastic, m::Model)

Update the values of a stochastic node according to its relationship with others in a model.

Returns the node with its values updated.

* `s` : stochastic node to update.

* `m` : model containing the node.
"""
function update!(s::AbstractStochastic, m::Model)
    #s.distr = s.eval(m)
    s = Stochastic(s, s.eval(m))
    s
end




#################### Distribution Methods ####################

function unlist(s::AbstractStochastic, transform::Bool = false)
    unlist(s, s.value, transform)
end

function unlist(s::AbstractStochastic, x::Real, transform::Bool = false)
    unlist(s, [x], transform)
end

function unlist(s::AbstractStochastic, x::AbstractArray, transform::Bool = false)
    transform ? unlist_sub(s.distr, link_sub(s.distr, x)) : unlist_sub(s.distr, x)
end

function relist(s::AbstractStochastic, x::AbstractArray, transform::Bool = false)
    relistlength(s, x, transform)[1]
end

function relistlength(s::AbstractVariate,
    x::AbstractArray, transform::Bool=false)
    value, n = relistlength_sub(s.distr, s, x)
    if transform
        u = invlink_sub(s.distr, value)
        return u, n
    else
        return value, n
    end
    #(transform ? invlink_sub(s.distr, value) : value, n)
end

function relistlength(s::TreeVariate, x::N, transform::Bool=false) where N<:GeneralNode
    value, n = relistlength_sub(s.distr, s, x)
    if transform
        return invlink_sub(s.distr, value)
    else
        return value, n
    end # if/else
end # relistlength

function logpdf(s::AbstractStochastic, transform::Bool=false)
    ll = logpdf(s, s.value, transform)
    
    s.lpdf = isa(ll, Float64) ? ll : ll.value
    ll
end

function logpdf(s::TreeVariate, transform::Bool=false)
    ll = logpdf(s, s.value, transform)
    s.lpdf = ll
    ll
end

function conditional_likelihood(s::AbstractStochastic, args...)
    conditional_likelihood(s, s.value, args...)
end

function conditional_likelihood(s::AbstractStochastic, x::AbstractArray, args...)
    logcond(s.distr, x, args...)
end


function pseudologpdf(s::AbstractStochastic, x::Union{Real, AbstractArray},
                      transform::Bool=false)
  logpdf(s, x, transform)
end

function rand(s::AbstractStochastic, x::Int64)
    rand(s.distr, x)
end

function gradlogpdf(s::Union{AbstractStochastic, AbstractLogical})
  gradlogpdf(s, s.value)
end

function gradlogpdf(s::TreeVariate, x::N, transform::Bool=false
                   ) where N <: GeneralNode
   gradlogpdf(s.distr, x)
end

function gradlogpdf(s::AbstractStochastic, x::AbstractArray)
    gradlogpdf_sub(s.distr, x)
end

function logpdf(s::AbstractStochastic, x::Union{AbstractArray, Real}, transform::Bool=false)
  logpdf_sub(s.distr, x, transform)
end

function logpdf(s::TreeVariate, x::GeneralNode, transform::Bool=false)
  logpdf_sub(s.distr, x, transform)
end

rand(s::AbstractStochastic) = rand_sub(s.distr, s.value)

#################### Core Model Functionality ####################

#################### Constructors ####################

"""
    Model(; iter::Integer=0, burnin::Integer=0,
      samplers::Vector{Sampler}=Sampler[], nodes...)

Construct a `Model` object that defines a model for MCMC simulation.

Returns a `Model` type object.

* `iter`: current iteration of the MCMC simulation.

* `burnin`: number of initial draws to be discarded as a burn-in sequence to allow for convergence.

* `samplers`: block-specific sampling functions.

* `nodes...`: arbitrary number of user-specified arguments defining logical and stochastic nodes in the model. Argument values must be `Logical` or `Stochastic` type objects. Their names in the model will be taken from the argument names.

"""
function Model(;
    iter::Integer = 0,
    burnin::Integer = 0,
    samplers::Vector{Sampler} = Sampler[],
    nodes...,
)

    nodedict = Dict{Symbol,AbstractDependent}()
    @inbounds for (key, value) in nodes
        isa(value, AbstractDependent) ||
            throw(ArgumentError("nodes are not all Dependent types"))
        node = deepcopy(value)
        isa(key, Symbol) || throw(ArgumentError("Something is wrong here"))
        node.symbol = key
        nodedict[key] = node
    end
    m = Model(nodedict, Sampler[], ModelState[], iter, burnin, false, false, -Inf64)
    dag = ModelGraph(m)
    dependentkeys = keys(m, :dependent)
    terminalkeys = keys(m, :stochastic)
    for v in vertices(dag.graph)
        vkey = dag.keys[v]
        if vkey in dependentkeys
            m[vkey].targets = intersect(dependentkeys, gettargets(dag, v, terminalkeys))
        end
    end
    setsamplers!(m, samplers)
end


#################### Indexing ####################

Base.getindex(m::Model, nodekey::Symbol) = m.nodes[nodekey]


function Base.setindex!(m::Model, value, nodekey::Symbol)
    node = m[nodekey]
    m.nodes[nodekey] = set_node(node, value)
end

function set_node(node::Logical, value::A)::Logical{A} where A
    Logical(node, value)
end

function set_node(node::Stochastic, value::A)::Stochastic{A} where A
    Stochastic(node, value)
end


function Base.setindex!(m::Model, values::Dict, nodekeys::Vector{Symbol})
    @inbounds for key in nodekeys
        m[key] = values[key]
    end
end

function Base.setindex!(m::Model, value, nodekeys::Vector{Symbol})
    length(nodekeys) == 1 || throw(BoundsError())
    m[first(nodekeys)] = value
end


Base.keys(m::Model) = collect(keys(m.nodes))
"""
    Base.keys(m::Model, ntype::Symbol, at...)

Extract the symbols (keys) for all existing nodes or for nodes of a specified type.

* `m` : model containing the nodes of interest.

* `ntype` : type of nodes to return. Options are
  * `:all` : all input, logical, and stochastic model nodes.

  * `:assigned` : nodes that have been assigned values.

  * `:block` : stochastic nodes being updated by the sampling block(s) `at::Integer=0` (default: all blocks).

  * `:dependent` : logical and stochastic (dependent) nodes in topologically sorted order.

  * `:independent` or `:input` : input (independent) nodes.

  * `:logical` : logical nodes.

  * `:monitor` : stochastic nodes being monitored in MCMC sampler output.

  * `:output` : stochastic nodes upon which no other stochastic nodes depend.

  * `:source` : nodes upon which the node `at::Symbol` or vector of nodes `at::Vector{Symbol}` depends.

  * `:stochastic` : stochastic nodes.

  * `:target` : topologically sorted nodes that depend on the sampling block(s) `at::Integer=0` (default: all blocks), node `at::Symbol` , or vector of nodes `at::Vector{Symbol}` .

* `at...` : additional positional arguments to be passed to the `ntype` options, as described above.
"""
function Base.keys(m::Model, ntype::Symbol, at...)
    ntype == :block ? keys_block(m, at...) :
    ntype == :all ? keys_all(m) :
    ntype == :assigned ? keys_assigned(m) :
    ntype == :dependent ? keys_dependent(m) :
    ntype == :independent ? keys_independent(m) :
    ntype == :input ? keys_independent(m) :
    ntype == :logical ? keys_logical(m) :
    ntype == :monitor ? keys_monitor(m) :
    ntype == :output ? keys_output(m) :
    ntype == :source ? keys_source(m, at...) :
    ntype == :stochastic ? keys_stochastic(m) :
    ntype == :target ? keys_target(m, at...) :
    throw(ArgumentError("unsupported node type $ntype"))
end

function keys_all(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for key in keys(m)
        node = m[key]
        if isa(node, AbstractDependent)
            push!(values, key)
            append!(values, node.sources)
        end
    end
    unique(values)
end

function keys_assigned(m::Model)::Array{Symbol}
    if m.hasinits
        values = keys(m)
    else
        values = Symbol[]
        @inbounds for key in keys(m)
            if !isa(m[key], AbstractDependent)
                push!(values, key)
            end
        end
    end
    values
end

function keys_block(m::Model, block::Integer = 0)::Array{Symbol}
    block == 0 ? keys_block0(m) : m.samplers[block].params
end

function keys_block0(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for sampler in m.samplers
        append!(values, sampler.params)
    end
    unique(values)
end

function keys_dependent(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for key in keys(m)
        if isa(m[key], AbstractDependent)
            push!(values, key)
        end
    end
    intersect(tsort(m), values)
end

function keys_independent(m::Model)::Array{Symbol}
    deps = Symbol[]
    @inbounds for key in keys(m)
        if isa(m[key], AbstractDependent)
            push!(deps, key)
        end
    end
    setdiff(keys(m, :all), deps)
end

function keys_logical(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for key in keys(m)
        if isa(m[key], AbstractLogical) || isa(m[key], TreeLogical)
            push!(values, key)
        end
    end
    values
end

function keys_monitor(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for key in keys(m)
        node = m[key]
        if isa(node, AbstractDependent) && !isempty(node.monitor)
            push!(values, key)
        end
    end
    values
end

function keys_output(m::Model)::Array{Symbol}
    values = Symbol[]
    dag = ModelGraph(m)
    @inbounds for v in vertices(dag.graph)
        vkey = dag.keys[v]
        if isa(m[vkey], AbstractStochastic) && !any_stochastic(dag, v, m)
            push!(values, vkey)
        end
    end
    values
end

keys_source(m::Model, nodekey::Symbol)::Array{Symbol} = m[nodekey].sources

function keys_source(m::Model, nodekeys::Vector{Symbol})::Array{Symbol}
    values = Symbol[]
    @inbounds for key in nodekeys
        append!(values, m[key].sources)
    end
    unique(values)
end

function keys_stochastic(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for key in keys(m)
        if isa(m[key], Stochastic)# || isa(m[key], TreeStochastic)
            push!(values, key)
        end
    end
    values
end

function keys_target(m::Model, block::Integer = 0)::Array{Symbol}
    block == 0 ? keys_target0(m) : m.samplers[block].targets
end

function keys_target0(m::Model)::Array{Symbol}
    values = Symbol[]
    @inbounds for sampler in m.samplers
        append!(values, sampler.targets)
    end
    intersect(keys(m, :dependent), values)
end

keys_target(m::Model, nodekey::Symbol)::Array{Symbol} = m[nodekey].targets

function keys_target(m::Model, nodekeys::Vector{Symbol})::Array{Symbol}
    values = Symbol[]
    @inbounds for key in nodekeys
        append!(values, m[key].targets)
    end
    intersect(keys(m, :dependent), values)
end


#################### Display ####################
"""
    Base.show(io::IO, m::Model)

Write a text representation of the model, nodes, and attributes to the current output stream.
"""
function Base.show(io::IO, m::Model)
    showf(io, m, Base.show)
end
"""
    showall(io::IO, m::Model)

Write a verbose text representation of the model, nodes, and attributes to the current output stream.
"""
function showall(io::IO, m::Model)
    showf(io, m, Base.showall)
end

function showf(io::IO, m::Model, f::Function)
    print(io, "Object of type \"$(summary(m))\"\n")
    width = displaysize()[2] - 1
    @inbounds for node in keys(m)
        print(io, string("-"^width, "\n", node, ":\n"))
        f(io, m[node])
        println(io)
    end
end


#################### Auxiliary Functions ####################

function names(m::Model, monitoronly::Bool)
    values = AbstractString[]
    @inbounds for key in keys(m, :dependent)
        if monitoronly
            if !isempty(m[key].monitor)
                nodenames = names(m, key)
                append!(values, nodenames)
            end
        else
            nodenames = names(m, key)
            append!(values, nodenames)
        end
    end
    values
end

function names(m::Model, nodekey::Symbol)
    node = m[nodekey]

    unlist(node, names(node))
end

function names(m::Model, nodekeys::Vector{Symbol})
    values = AbstractString[]
    @inbounds for key in nodekeys
        append!(values, names(m, key))
    end
    values
end


function update_likelihood!(m::Model, likelihood::Float64)
    m.likelihood = likelihood
    m
end

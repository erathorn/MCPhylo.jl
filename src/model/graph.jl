#################### Model Graph ####################

function ModelGraph(m::Model)
    allkeys = keys(m, :all)
    g = DiGraph(length(allkeys))
    lookup = Dict(allkeys[i] => i for i = 1:length(allkeys))
    for key in keys(m)
        node = m[key]
        if isa(node, AbstractDependent)
            for src in node.sources
                add_edge!(g, lookup[src], lookup[key])
            end
        end
    end
    ModelGraph(g, allkeys)
end


#################### Display ####################
"""
    draw(m::Model; filename::AbstractString="")

Draw a GraphViz DOT-formatted graph representation of model nodes and their relationships.

The model drawn to an external file or standard output. Stochastic, logical, and input nodes will be represented by ellipses, diamonds, and rectangles, respectively. Nodes that are unmonitored in MCMC simulations will be gray-colored.

* `m` : model for which to construct a graph.

* `filename` : external file to which to save the resulting graph, or an empty string to draw to standard output (default). If a supplied external file name does not include a dot (`.`), the file extension `.dot` will be appended automatically.
"""
function draw(m::Model; filename::AbstractString = "")
    dot = graph2dot(m)
    if length(filename) == 0
        print(dot)
    else
        if something(findfirst(isequal('.'), filename), 0) == 0
            filename = string(filename, ".dot")
        end
        f = open(filename, "w")
        write(f, dot)
        close(f)
    end
end

graph(m::Model) = ModelGraph(m)
"""
    graph2dot(m::Model)

Draw a GraphViz DOT-formatted graph representation of model nodes and their relationships.

A character string representation of the graph suitable for in-line processing. Stochastic, logical, and input nodes will be represented by ellipses, diamonds, and rectangles, respectively. Nodes that are unmonitored in MCMC simulations will be gray-colored.

* `m` : model for which to construct a graph.
"""
function graph2dot(m::Model)
    dag = ModelGraph(m)
    io = IOBuffer()
    write(io, "digraph MCPhyloModel {\n")
    deps = keys(m, :dependent)
    for v in vertices(dag.graph)
        attr = Tuple{AbstractString,AbstractString}[]
        vkey = dag.keys[v]
        if vkey in deps
            node = m[vkey]
            if isa(node, AbstractLogical)
                push!(attr, ("shape", "diamond"))
            elseif isa(node, AbstractStochastic)
                push!(attr, ("shape", "ellipse"))
            end
            if isempty(node.monitor)
                push!(attr, ("style", "filled"), ("fillcolor", "gray85"))
            end
        else
            push!(attr, ("shape", "box"), ("style", "filled"), ("fillcolor", "gray85"))
        end
        write(io, "\t\"")
        write(io, vkey)
        write(io, "\" [")
        write(io, join(map(x -> "$(x[1])=\"$(x[2])\"", attr), ", "))
        write(io, "];\n")
        for t in outneighbors(dag.graph, v)
            write(io, "\t\t\"")
            write(io, vkey)
            write(io, "\" -> \"")
            write(io, dag.keys[t])
            write(io, "\";\n")
        end
    end
    write(io, "}\n")
    String(take!(copy(io)))
end


#################### Auxiliary Functions ####################

function any_stochastic(dag::ModelGraph, v::Int, m::Model)
    found = false
    for t in outneighbors(dag.graph, v)
        tkey = dag.keys[t]
        if isa(m[tkey], AbstractStochastic) || any_stochastic(dag, t, m)
            found = true
            break
        end
    end
    found
end

function gettargets(dag::ModelGraph, v::Int, terminalkeys::Vector{Symbol})
    values = Symbol[]
    for t in outneighbors(dag.graph, v)
        tkey = dag.keys[t]
        push!(values, tkey)
        if !(tkey in terminalkeys)
            values = union(values, gettargets(dag, t, terminalkeys))
        end
    end
    values
end

function tsort(m::Model)
    dag = ModelGraph(m)
    dag.keys[topological_sort_by_dfs(dag.graph)]
end

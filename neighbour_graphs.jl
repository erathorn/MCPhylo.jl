using Distances, LinearAlgebra, LightGraphs, MetaGraphs

function create_neighbour_list(dm::Array{Float64,2}, threshold::Int64, languages::Vector{String})
    neighbours = Vector()
    for (i, l1) in enumerate(languages)
        for (j, l2) in enumerate(languages)
            if (dm[i,j] <= threshold) && (dm[i,j] != 0)
                push!(neighbours, [l1,l2])
            end
        end
    end
    return neighbours
end

function create_nmat(dm::Array{Float64,2}, threshold::Int64)
    n, m = size(dm)
    nmat = zeros(n, m)
    idx_list = Vector(1:n)
    for (i, l1) in enumerate(idx_list)
        for (j, l2) in enumerate(idx_list)
            if (dm[i,j] <= threshold) && (dm[i,j] != 0)
                nmat[i,j] = 1
            end
        end
    end
    return nmat
end

function replace_diagonal!(nmat::Matrix)
    for i in diagind(nmat, 0)
        nmat[i] == 1 ? nmat[i] = 0 : nothing
    end
    return nmat
end

function create_linguistic_nmat(d::DataFrame)
    nmat = [Int(g1==g2) for g1 in d.Genus, g2 in d.Genus]
    replace_diagonal!(nmat)
    return nmat
end

# to do
function create_weighted_graph(x)
    return x
end


"""
To set the value of y for each vertex in the graph, convert
the SimpleGraph (ng) to a MetaGraph (mg).
Use set_prop! to set a property for a vertex or edge;
set_yvals! does t[his for the whole graph.
"""

function set_yvals!(g::MetaGraph, Y::Vector) # or: rewrite to take df as input
    for (i, y) in enumerate(Y)
        set_prop!(g, i, :y, y)
    end
    return g
end

"""
The functions below find neighbouring pairs for vertex v
that are concordant for trait y.
"""

function concordant_sum(X::Array{Int64,1}, y::Symbol, g::MetaGraph, v::Int64)
    concordant_sum = 0
    for n in neighbors(g, v)
#        if has_prop(mg, :x, v) - build this in to throw an error if prop not defined for g, v
        if get_prop(g, n, :y) == get_prop(g, v, :y)
            concordant_pairs += 1 # or another way to measure the pairs?
        end
    end
    return concordant_sum
end

# Example

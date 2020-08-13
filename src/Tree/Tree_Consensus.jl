"""
    are_compatible(cluster::Vector{T}, tree::T) where T<:AbstractNode

Check if a cluster is compatible with a tree
"""
function are_compatible(cluster::Vector{T}}, tree::T) where T<:AbstractNode
    lca = find_lca(tree, cluster)
    children = lca.children
    for child in children
        leaves = get_leaves(child)
        if !(length(intersect(leaves, cluster)) in (0, length(leaves)))
            return false
        end #if
    end # end
    return true
end

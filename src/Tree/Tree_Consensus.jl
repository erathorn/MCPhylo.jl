"""
    are_compatible(cluster_names::Vector{String}, tree::T) where T<:AbstractNode

Check if a cluster of nodes, derived from the node names is compatible with a tree
"""
function are_compatible(cluster::Vector{String}, tree::T)::Bool where T<:AbstractNode
    tree_leaves = get_leaves(tree)
    leaf_names = [leaf.name for leaf in tree_leaves]
    if length(intersect(leaf_names, cluster)) != length(cluster)
        throw(ArgumentError("The cluster contains non-leaf nodes."))
    end # if
    lca = find_lca(tree, cluster)
    children = lca.children
    for child in children
        leaves = get_leaves(child)
        leaves = [leaf.name for leaf in leaves]
        if !(length(intersect(leaves, cluster)) in (0, length(leaves)))
            return false
        end #if
    end # end
    return true
end

"""
    are_compatible(cluster::Vector{T}, tree::T) where T<:AbstractNode

Check if a cluster of nodes is compatible with a tree
"""
function are_compatible(cluster::Vector{T}, tree::T)::Bool where T<:AbstractNode
    cluster = [node.name for node in cluster]
    are_compatible(cluster, tree)
end

function days_algorithm(cluster1::T, cluster2::T) where T<:AbstractNode

end

"""
    are_compatible(cluster_names::Vector{String}, tree::T) where T<:AbstractNode

Check if a cluster of nodes, derived from the node names is compatible with a tree
"""
function are_compatible(cluster::Vector{String}, tree::T)::Bool where T<:AbstractNode
    leaf_names = [leaf.name for leaf in get_leaves(tree)]
    if length(intersect(leaf_names, cluster)) != length(cluster)
        throw(ArgumentError("The cluster contains non-leaf nodes."))
    end # if
    lca = find_lca(tree, cluster)
    children = lca.children
    for child in children
        leaves = [leaf.name for leaf in get_leaves(child)]
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


"""
    majority_consensus_tree(trees::Vector{T}, tree::T)::T where T<:AbstractNode

Construct the majority rule consensus tree from a set of trees
"""
function majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode
    all_leaves = [get_leaves(root) for root in trees]
    leaf_names = sort([leaf.name for leaf in all_leaves[1]])
    for leaves in all_leaves[2:end]
        if leaf_names != sort([leaf.name for leaf in leaves])
            throw(ArgumentError("Input Trees need to have identical leaf names"))
        end
    end
    # TODO: need Day's Algorithm next
    return trees[1] # to make it not crash
end

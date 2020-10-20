"""
    prune_tree!(root::T, node_names::Vector{String})::Nothing

This function removes specific nodes including their descendants from the tree.
"""
function prune_tree!(root::T, node_names::Vector{String})::Nothing where T<:AbstractNode
    if root.name in node_names
        throw(ArgumentError("trying to prune root, please set root to nothing instead"))
    else
        for child in root.children
            if child.name in node_names
                remove_child!(root, child)
            end
        end
        # recursively check all nodes of the tree
        for child in root.children
            prune_tree!(child, node_names)
        end
    end

end


"""
    prune_tree(root::T, node_names::Vector{T})::T where T<:AbstractNode

This function returns a copy of a tree with specific nodes including their
descendants removed
"""
function prune_tree(root::T, nodes::Vector{T})::T where T<:AbstractNode
    # copy the tree and call the inplace version of the function on the copy
    if root.name in node_names
        throw(ArgumentError("trying to prune root, please set root to nothing instead"))
    end
    copyroot = deepcopy(root)
    prune_tree!(copyroot, nodes)
    return copyroot
end


"""
    prune_tree!(root::T, node_names::Vector{String})::Nothing where T<:AbstractNode

In-place version of prune_tree
"""
function prune_tree!(root::T, node_names::Vector{String})::Nothing where T<:AbstractNode
    nodes = post_order(root)
    names = [node.name for node in nodes]
    nodes_to_prune = Vector{Node}()
    for name in node_names
        indeces = findall(x->x == name, names)
        if length(indeces) > 1
            throw(ArgumentError("Multiple nodes are named \"$name\""))
        elseif length(indeces) == 0
            print("Warning: No node named \"$name\" in tree")
        else
            push!(nodes_to_prune, nodes[indeces[1]])
        end # if/else
    end # for
    if length(nodes_to_prune) == 0
        throw(ArgumentError("None of the node names correspond to a node in the tree"))
    end # if
    prune_tree!(root, nodes_to_prune)
end


"""
    prune_tree!(root::T, node_names::Vector{T})::Nothing where T<:AbstractNode

In-place version of prune_tree
"""
function prune_tree!(root::T, nodes::Vector{T})::Nothing where T<:AbstractNode
    if root in nodes
        throw(ArgumentError("trying to prune root, please set root to nothing instead"))
    else
        for node in nodes
            remove_child!(node.mother, node)
        end
    end
end

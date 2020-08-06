"""
    prune_tree!(root::T, node_names::Vector{String})

This function removes specific nodes and their daughters from a tree
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
    prune_tree(root::T, node_names::Vector{String})::T where T<:AbstractNode

This function returns a copy of a tree with specific nodes including their
descendants removed
"""
function prune_tree(root::T, node_names::Vector{String})::T where T<:AbstractNode
    # copy the tree and call the inplace version of the function on the copy
    if root.name in node_names
        throw(ArgumentError("trying to prune root, please set root to nothing instead"))
    end
    copyroot = deepcopy(root)
    prune_tree!(copyroot, node_names)
    return copyroot
end

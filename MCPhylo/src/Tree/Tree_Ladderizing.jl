function ladderize_tree!(root::T,ascending::Bool=true)::Nothing where T<:AbstractNode
    if root.nchild == 0
        return
    end
    ndescendants = Array{Int64,1}(undef, length(root.children))
    for (index, child) in enumerate(root.children)
        ndescendants[index] = length(post_order(child))
    end
    perm = sortperm(ndescendants)
    if ascending == false
        reverse!(perm)
    end
    root.children = root.children[perm]

    for child in root.children
        ladderize_tree!(child, ascending)
    end
end

function ladderize_tree(root::T, ascending::Bool=true)::T where T<:AbstractNode
    copyroot = deepcopy(root)
    ladderize_tree!(copyroot)
    return copyroot
end

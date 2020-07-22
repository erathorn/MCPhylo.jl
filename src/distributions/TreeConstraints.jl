
function topological(tree::N, constraints::Dict) where N<:GeneralNode
    for key in keys(constraints)
        lca = find_lca(tree, constraints[key])
        lca.root && return false
    end
    true
end


# topological constraints fallback
function topological(tree::N, constraints::Missing) where N<:GeneralNode
    true
end

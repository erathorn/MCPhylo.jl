
function topological(tree::Node, constraints::Dict)
    for key in keys(constraints)
        lca = find_lca(tree, constraints[key])
        lca.root && return false
    end
    true
end


# topological constraints fallback
function topological(tree::Node, constraints::Missing)
    true
end

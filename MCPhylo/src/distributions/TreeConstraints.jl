
function topological(tree::Node, constraints::Dict)
    for key in keys(constraints)
        println(key)
    end
end


# topological constraints fallback
function topological(tree::Node, constraints::Missing)
    true
end

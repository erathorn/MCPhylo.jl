
#################### Post order traversal ####################
"""
    post_order(root::T, traversal::Vector{T})::Vector{T} where T<:GeneralNode

This function performs a post order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the post order traversal.

Returns vector of Nodes.

* `root` : root Node of tree.

* `traversal` : vector of Nodes; nodes are pushed to this vector as the tree is traversed.
"""
function post_order(root::T, traversal::Vector{T}) where T<:GeneralNode
   if root.nchild != 0
        for child in root.children
            post_order(child, traversal)
        end
   end # if
   push!(traversal, root)
   return traversal
end # function post_order_trav


"""
    post_order(root::T)::Vector{T} where T<:GeneralNode

This function does post order traversal. Only the root node needs to be supplied.

Returns vector of Nodes.

* `root` : root Node of tree.
"""
function post_order(root::T)::Vector{T} where T<:GeneralNode
    t::Vector{T} = []
    post_order(root, t)
    return t
end # function post_order

"""
    get_leaves(root::T, traversal::Vector{T})::Vector{T} where T<:GeneralNode

This function returns leaf nodes of a tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used.

Returns a vector of leaf Nodes.

* `root` : root Node of tree to traverse.

* `traversal` : vector of Nodes; leaf Nodes, once found, are pushed to this vector.

"""
function get_leaves(root::T, traversal::Vector{T})::Vector{T} where T<:GeneralNode
   if root.nchild != 0
        for child in root.children
            get_leaves(child, traversal)
        end
   else
       push!(traversal, root)
   end # if

   return traversal
end # function post_order_trav


"""
    get_leaves(root::T)::Vector{T} where T<:GeneralNode

This function returns the leaves of a tree. Only the root node needs to be supplied.

Returns vector of leaf Nodes.

* `root` : root Node of tree.
"""
function get_leaves(root::T)::Vector{T} where T<:GeneralNode
    t::Vector{T} = []
    get_leaves(root, t)
    return t
end # function post_order




#################### Pre order traversal ####################

"""
    pre_order(root::T, traversal::Vector{T})::Vector{T} where T<:GeneralNode

This function performs a pre order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the pre order traversal.

Returns vector of Nodes.

* `root` : root Node of tree.

* `traversal` : vector of Nodes; nodes are pushed to this vector as the tree is traversed.
"""
function pre_order(root::T, traversal::Vector{T})::Vector{T} where T<:GeneralNode
    push!(traversal, root)
    if root.nchild != 0
        for child in root.children
            pre_order(child, traversal)
        end
    end # if
    return traversal
end # function pre_order!


"""
    pre_order(root::T)::Vector{T} where T<:GeneralNode

This function does pre order traversal. Only the root node needs to be supplied.

Returns vector of Nodes.

* `root` : root Node of tree.
"""
function pre_order(root::T)::Vector{T} where T<:GeneralNode
    t::Vector{T} = []
    pre_order(root, t)
    return t
end # function pre_order

#################### Level order traversal ####################

"""
    level_order(node::T)::Array{T} where T<:GeneralNode

This function does level order traversal. Only the root node needs to be supplied.

Returns Array of Nodes.

* `node` : root Node of tree.
"""
function level_order(node::T)::Array{T} where T<:GeneralNode
    level = 1
    stack::Array{T} = []
    while level_traverse(node, level, stack)
        level += 1
    end # while
    stack
end # function level_order

"""
    level_traverse(node::T, level::Int64, stack::Array{T})::Bool where T <:GeneralNode

This function traverses a level of the tree specified through `node`. The level
is specified via the `level` argument and the nodes visited are stored in the
`stack`.
This function is intended as the internal worker for the level_order function.
"""
function level_traverse(node::T, level::Int64, stack::Array{T})::Bool where T <:GeneralNode

    if level == 1
        # level which needs to be traversed right now
        push!(stack, node)
        return true
    else
        # move down the tree to the correct level
        boolqueue = [false] # this is used to look for the correct level
        for child in node.children
            push!(boolqueue, level_traverse(child, level-1, stack))
        end # for
        return reduce(|, boolqueue)
    end # if
end # function level_traverse

#=
my_tree:
- Julia version: 1.3.1
- Author: erathorn
- Date: 2019-05-07
=#

#TODO: Automate export of automatically genereated funtions

"""
    add_child!(mother_node::Node, child::Node)

This function adds a child to the mother node.
The arity of the mother node is increased by `1` and the root
status of the child is set to `False`.
"""
function add_child!(mother_node::Node, child::Node)
    push!(mother_node.children, child)
    child.mother = mother_node
    mother_node.nchild += 1
    child.root = false
    mother_node.initialized=true
end # function add_child


"""
    remove_child!(mother_node::Node, index::int)Node

This function removes a child from the list of nodes which are daughters of this
node. The removed node is returned.
"""
function remove_child!(mother_node::Node, left::Bool)::Node
    @assert (mother_node.nchild === 2)
    if left
        rv = popfirst!(mother_node.children)
        rv.mother = missing
    else
        rv = pop!(mother_node.children)
        rv.mother = missing
    end # end if
    mother_node.nchild -= 1
    return rv
end # function

"""
    remove_child!(mother_node::Node, index::int)Node

This function removes a child from the list of nodes which are daughters of this
node. The removed node is returned.
"""
function remove_child!(mother_node::Node, child::Node)::Node
    ind = findfirst(x->x==child, mother_node.children)
    deleteat!(mother_node.children, ind)
    child.mother = missing

    mother_node.nchild -= 1
    return child
end # function


function tree_from_leaves(leaf_nodes::Vector{String},node_size::Int, final_length::Int64)::Tuple{Vector{Node}, Int}
    my_node_list::Array{Node,1} = []

    # first create a list of leaf nodes
    for node_name in leaf_nodes
        nn = Node(node_name, data=zeros(Float64, (2, node_size)))

        push!(my_node_list,nn)
    end # for

    # Internal nodes are created using integers as names.
    temp_name::Int = length(my_node_list)+1

    # shuffle the node list to get a random tree
    Random.shuffle!(my_node_list)

    while length(my_node_list) > final_length
        # get two nodes
        # create a new mother node to which the two first nodes are added as children
        # add the new mother node to the list and reshuffle
        first_child::Node = pop!(my_node_list)
        first_child.inc_length = rand(Uniform(0.0015,1))#*0.1
        second_child::Node = pop!(my_node_list)
        second_child.inc_length = rand(Uniform(0.0015,1))
        curr_node::Node = Node(string(temp_name), data=zeros(Float64, (2, node_size)))

        add_child!(curr_node, first_child)
        add_child!(curr_node, second_child)
        push!(my_node_list, curr_node)
        temp_name += 1
        Random.shuffle!(my_node_list)
    end # while

    return my_node_list, temp_name
end


"""
    create_tree_from_leaves(leaf_nodes::Vector{T})::Node

This function creates a  random binary tree from a list of leaf nodes.
The root node as access point for the tree is returned.
"""
function create_tree_from_leaves_bin(leaf_nodes::Vector{String}, node_size::Int)::Node

    my_node_list, temp_name = tree_from_leaves(leaf_nodes, node_size ,2)

    root::Node = Node(string(temp_name), data=zeros(Float64, (2, node_size)))

    lchild = pop!(my_node_list)
    lchild.inc_length = rand(Uniform(0.0015,1))

    rchild = pop!(my_node_list)
    rchild.inc_length = rand(Uniform(0.0015,1))
    add_child!(root, lchild)
    add_child!(root, rchild)


    set_binary!(root)
    number_nodes!(root)

    return root
end # function create_tree_from_leaves



"""
    create_tree_from_leaves(leaf_nodes::Vector{T})::Node

This function creates a  random binary tree from a list of leaf nodes.
The root node as access point for the tree is returned.
"""
function create_tree_from_leaves(leaf_nodes::Vector{String}, node_size::Int64 = 1; cu::Bool=false)::N where N <: AbstractNode

    my_node_list, temp_name = tree_from_leaves(leaf_nodes, 3)

    root::Node = Node(string(temp_name), data=zeros(Float64, (2, node_size)))
    lchild = pop!(my_node_list)
    lchild.inc_length = rand(Uniform(0.0015,1))
    mchild = pop!(my_node_list)
    mchild.inc_length = rand(Uniform(0.0015,1))
    rchild = pop!(my_node_list)
    rchild.inc_length = rand(Uniform(0.0015,1))
    add_child!(root, lchild, true)
    add_child!(root, rchild, false)
    add_child!(root, mchild, false, true)

    set_binary!(root)
    number_nodes!(root)

    return root
end # function create_tree_from_leaves


function create_tree_from_leaves_cu(leaf_nodes::Vector{String}, node_size::Int64 = 1)::Node_cu
    my_node_list::Array{Node_cu,1} = []

    # first create a list of leaf nodes
    for node_name in leaf_nodes
        nn =  Node_cu(node_name, zeros(Float64, (2, node_size)),missing, missing, missing, 0, true, 0.0, "0", 0, 0.0)
        push!(my_node_list,nn)
    end # for

    # Internal nodes are created using integers as names.
    temp_name::Int = length(my_node_list)+1

    # shuffle the node list to get a random tree
    Random.shuffle!(my_node_list)

    while length(my_node_list) != 1
        # get two nodes
        # create a new mother node to which the two first nodes are added as children
        # add the new mother node to the list and reshuffle
        first_child::Node = pop!(my_node_list)
        first_child.inc_length = rand(Uniform(0,1))
        second_child::Node = pop!(my_node_list)
        second_child.inc_length = rand(Uniform(0,1))
        curr_node::Node = Node_cu(string(temp_name), zeros(Float64, (2, node_size)), missing, missing, missing, 0, true, 0.0, "0", 0,0.0)
        add_child!(curr_node, first_child, true)
        add_child!(curr_node, second_child, false)
        push!(my_node_list, curr_node)
        temp_name += 1
        Random.shuffle!(my_node_list)
    end # while
    root = pop!(my_node_list)
    set_binary!(root)
    number_nodes!(root)

    return root
end # function create_tree_from_leaves



# legacy wrapper
function rescale_length(root::T) where T<:AbstractNode
    force_ultrametric(root)
end


"""
    force_ultrametric!(root::T) where T<:AbstractNode

Force an ultrametric version of the tree.
"""
function force_ultrametric!(root::T) where T<:AbstractNode
    po::Vector{T} = post_order(root)
    node2max_depth = zeros(UInt32, length(po))
    for node in po
        if node.nchild != 0
            mv = -1
            for child in node.children
                if node2max_depth[child.num] > mv
                    mv = node2max_depth[child.num]
                end # if
            end # for
            node2max_depth[node.num] = mv+1
        else
            node2max_depth[node.num] = 1
        end # if
    end # for

    node2dist = zeros(Float64, size(node2max_depth))
    tl = tree_height(root)
    nblv = zeros(Float64,length(po))
    for node in level_order(root)
        if node.root != true
            nv = (tl - node2dist[node.mother.num])/node2max_depth[node.num]

            nblv[node.num] = nv
            node2dist[node.num] = nv + node2dist[node.mother.num]
        end # if
    end # for

    set_branchlength_vector!(root, nblv)
end # function force_ultrametric!


#################### Tree length & height ####################

"""
    tree_length(root::Node)::Float64

This function calculates the tree_length.
"""
function tree_length(root::T)::Float64  where T<:AbstractNode
    return tree_length(root, 0.0)
end # function tree_length

"""
    tree_length(root::T, tl::Float64)::Float64 where T<:AbstractNode

This function does the internal tree length recursion
"""
function tree_length(root::T, tl::Float64)::Float64 where T<:AbstractNode

    if length(root.children) != 0
        for child in root.children
            tl = tree_length(child, tl)
        end
    end # if
    if root.root !== true
        tl += root.inc_length
    end

    tl
end # function tree_length



"""
    tree_height(root::Node)::Float64

This function calculates the tree height.
"""
function tree_height(root::T)::Float64  where T<:AbstractNode
    node_height(root)
    return root.height
end

"""
    node_height(root::T, mv::Float64)::Float64  where T<:AbstractNode

Calculate the height of a node.
"""
function node_height(root::T)  where T<:AbstractNode

    if root.nchild != 0
        for node in root.children
            node_height(node)
        end
        root.height = maximum([child.inc_length+child.height for child in root.children])
    else
        root.height = 0.0
    end
end # function node_height

function node_height_vec(root::T, vec::Vector{N})  where {T<:AbstractNode, N<:Real}

    if root.nchild != 0
        for node in root.children
            node_height_vec(node, vec)
        end
        root.height = maximum([child.inc_length+child.height for child in root.children])
    else
        root.height = 0.0
    end
    vec[root.num] = root.height
end # function node_height


function node_height_vec(root::T)::Vector{Float64} where T<:AbstractNode
    t = zeros(length(post_order(root)))
    node_height_vec(root, t)
    t
end # function node_height


function node_distance(tree::T, node1::T, node2::T)::Float64 where T<:AbstractNode
    lca = find_lca(tree, node1, node2)
    path_length(lca, node1)+path_length(lca,node2)
end

function get_path(ancestor::T, descendant::T)::Vector{Int64} where T<:AbstractNode
    path::Vector{Int64} = []
    while descendant.num != ancestor.num
        push!(path, descendant.num)
        descendant = descendant.mother
    end
    path
end


"""
    path_length(ancestor::Node, descendant::Node)::Float64

Note: The function assumes there is an ancestral relationship between the two nodes.

This function calculates the length of the path separating the ancestor from the
offspring node. The function follows the path specified through the binary
description of the node.
"""
function path_length(ancestor::T, descendant::T)::Float64  where T<:AbstractNode
    l::Float64 = 0.0

    while descendant != ancestor
        l += descendant.inc_length
        descendant = descendant.mother
    end # while
    return l
end #function path_length


"""
    get_sister(root::Node, node::Node)::Node

This function gets the sister of `node`. It does so by looking for the respective
binary representation of the sister.
"""
@inline function get_sister(node::T)::T  where T<:AbstractNode
    node.mother.children[findfirst(y-> y!=node, node.mother.children)]
end # function


"""
    get_mother(root::Node, node::Node)::Node

This function gets the mother of `node`. It does so by looking for the respective
binary representation of the mother node.
"""
@inline function get_mother(node::T)::T  where T<:AbstractNode
    return node.mother
end # function

"""
    set_binary!(root::Node)

Assign a binary representation to each node, which specifies the path from the
root to this node via the binary representation of the node.
A left turn is a 1 in binary and a right turn a 0.
"""
function set_binary!(root::T)  where T<:AbstractNode
    if root.root
        root.binary = "1"
    end # if
    if root.nchild != 0
        for (ind, node) in enumerate(root.children)
            ind -= 1
            node.binary = string(root.binary, ind)
            set_binary!(node)
        end

    end # if
end # function set_binary

"""
    number_nodes!(root::Node)::Nothing

This function assigns a unique, sequential number to each node.
"""
function number_nodes!(root::T)::Nothing  where T<:AbstractNode
    for (index, value) in enumerate(post_order(root))
        value.num = index
    end # for
end # fuction number_nodes


"""
    random_node(root::Node)::Node

This function returns a random node from the tree.
"""
function random_node(root::T)::T  where T<:AbstractNode
    post_order_trav = post_order(root)
    return rand(post_order_trav)
end # function random_node



#################### Vector of branch lengths: get & set ####################

"""
    get_branchlength_vector(root::N)::Vector{T}  where {N <:AbstractNode, T<:Real}

Get the vector of branch lengths of the tree.
"""
function get_branchlength_vector(root::N)::Vector{Float64}  where {N <:AbstractNode}
    get_branchlength_vector(root, root.blv)
    return root.blv
end # function get_branchlength_vector

"""
    get_branchlength_vector(t::TreeStochastic)

Get the vector of branch lengths of the tree.
"""
function get_branchlength_vector(t::TreeStochastic)::Vector{Float64}
    get_branchlength_vector(t.value)
end # function

"""
    get_branchlength_vector(root::N, out_vec::Vector{T}) where {N<:AbstractNode, T<:Real}

Do post order traversal to retrieve a vector of branch lengths.
"""
function get_branchlength_vector(root::N, out_vec::Vector{T})::Nothing where {N<:Node{<:Real,<:Real,<:Real,<:Integer}, T<:Real}
    for child in root.children
        get_branchlength_vector(child, out_vec)
    end
    if !root.root
        out_vec[root.num] = root.inc_length
    end
    nothing
end

"""
    get_branchlength_vector(root::N, out_vec::Vector{T}) where {N<:AbstractNode, T<:Real}

Do post order traversal to retrieve a vector of branch lengths.
"""
function get_branchlength_vector(root::N, vec::Nothing)::Vector{Float64} where {N<:AbstractNode, T<:Real}
    root.blv = Vector{Float64}(undef, length(post_order(root))-1)
    vec = root.blv
    get_branchlength_vector(root, vec)
    return vec
end



"""
    set_branchlength_vector!(t::TreeStochastic, blenvec::Array{T}) where T <: Real

Get the vector of branch lengths of the tree.
"""
function set_branchlength_vector!(t::TreeStochastic, blenvec::Array{T}) where T <: Real
    set_branchlength_vector!(t.value, blenvec)
end # function

"""
    set_branchlength_vector!(t::TreeStochastic, blenvec::Array{T}) where T <: Real

Get the vector of branch lengths of the tree.
"""
function set_branchlength_vector!(t::TreeStochastic, blenvec::ArrayStochastic)
    set_branchlength_vector!(t.value, blenvec.value)
end # function

"""
    set_branchlength_vector!(root::Node, blenvec::Vector{Float64})::Node

This function sets the branch lengths of a tree to the values specified in blenvec.
"""
function set_branchlength_vector!(root::N, blenvec::Array{T}) where {N<:AbstractNode, T<:Real}
    any(0 .> blenvec) && throw("this should never happen")
    for child in root.children
        set_branchlength_vector!(child, blenvec)
    end

    @views if root.root !== true
        root.inc_length = blenvec[root.num]
    end
    nothing
end # function set_branchlength_vector!


"""
    get_sum_seperate_length!(root::Node)::Vector{Float64}

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(root::T)::Vector{Float64}  where T<:AbstractNode
    return get_sum_seperate_length!(post_order(root))
end # function get_sum_seperate_length!


"""
    get_sum_seperate_length!(root::Node)::Vector{Float64}

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(post_order::Vector{T})::Vector{Float64}  where T<:AbstractNode
    res_int::Float64 = 0.0
    res_leave::Float64 = 0.0
    res_int_log::Float64 = 0.0
    res_leave_log::Float64 = 0.0
    @simd for node in post_order
        if node.nchild !== 0
            # internal branches
            if !node.root
                res_int += node.inc_length
                res_int_log += log(node.inc_length)
            end
        else
            # branches leading to leaves
            res_leave += node.inc_length
            res_leave_log += log(node.inc_length)
        end # if
    end # for
    return [res_int, res_leave, res_int_log, res_leave_log]
end # function get_sum_seperate_length!

function internal_external_map(t::TreeStochastic)::Vector{Int64}
    internal_external_map(t.value)
end

function internal_external_map(root::T)::Vector{Int64}  where T<:AbstractNode
    internal_external_map(post_order(root))
end

function internal_external_map(post_order::Vector{T})::Vector{Int64}  where T<:AbstractNode
    my_map::Vector{Int64} = zeros(Int64, length(post_order)-1)
    for node in post_order
        if !node.root
            if node.nchild != 0
                my_map[node.num] = 1
            end
        end
    end
    return my_map
end

function internal_external(root::T)::Vector{Int64}  where T<:AbstractNode
    v = root.IntExtMap
    if v === nothing
        v = internal_external_map(root)
        root.IntExtMap = v
    end
    v
end


function find_lca(tree::T, node_l::Array{String, 1})::T  where T<:AbstractNode
    find_lca(tree, [find_by_name(tree, i) for i in node_l])
end

function find_lca(tree::T, node_l::Array{T})::T  where T<:AbstractNode
    @assert length(node_l) > 0
    if length(node_l) === 1
        return node_l[1]
    else
        n1 = popfirst!(node_l)
        n2 = popfirst!(node_l)
        lca = find_lca(tree, n1, n2)
        while length(node_l) !== 0
            n1 = popfirst!(node_l)
            lca = find_lca(tree, lca, n1)
        end
        return lca
    end
end

function find_lca(tree::T, node1::T, node2::T)::T  where T<:AbstractNode
    nb = lcp(node1.binary, node2.binary)
    find_binary(tree, nb)
end

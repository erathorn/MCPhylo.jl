#=
my_tree:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#



#module Tree_Basics

#using Markdown
#using Random
#using Distributions


#export create_tree_from_leaves, post_order, tree_length, tree_height,
#       path_length, get_leaves, Node, add_child!, set_binary!, remove_child!,
#       random_node, move!, find_by_name, find_by_binary, find_by_root, get_branchlength_vector!


#TODO: Automate export of automatically genereated funtions

"""
    Node

This data type holds the basic Node structure. The type T is used to specify the type of the data
stored in the node.

* If `nchild` is `0` the Node is a leaf node.
* If `root` is `False` the Node is a child of another node.
* `inc_length` specifies the length of the incomming branch.
* `binary` specifies the path from the root to the Node. `1` and `0` represent left and right turns respectively.
"""
mutable struct Node
    name::String
    data::Array{Float64,2}
    child::Vector{Node}
    nchild::Int64
    root::Bool
    inc_length::Float64
    binary::String
end # struct Node


Node() = Node(String[],Float64[], Node[], 0, true, 0.0, "0")


"""
    add_child(mother_node::Node, child::Node)

This function adds a child to the mother node.
The arity of the mother node is increased by `1` and the root
status of the child is set to `False`.
"""
function add_child!(mother_node::Node, child::Node)
    push!(mother_node.child, child)
    mother_node.nchild += 1
    child.root = false
end # function add_child

"""
    remove_child!(mother_node::Node, index::int)Node

This function removes a child from the list of nodes which are daughters of this
node. The removed node is returned.
"""
function remove_child!(mother_node::Node, index::Int)
    mother_node.nchild -= 1
    if index == 1
        return popfirst!(mother_node.child)
    else
        return pop!(mother_node.child)
    end # end if
end # function



"""
    create_tree_from_leaves(leaf_nodes::Vector{T})::Node

This function creates a  random binary tree from a list of leaf nodes.
The root node as access point for the tree is returned.
"""
function create_tree_from_leaves(leaf_nodes::Vector{String}, node_size::Int64 = 1)::Node
    my_node_list::Array{Node,1} = []

    # first create a list of leaf nodes
    for node_name in leaf_nodes
        push!(my_node_list, Node(node_name, zeros(Float64, (2, node_size)), Node[], 0, true, 0.0, "0"))
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
        curr_node::Node = Node(string(temp_name), zeros(Float64, (2, node_size)), Node[], 0, true, 0.0, "0")
        add_child!(curr_node, first_child)
        add_child!(curr_node, second_child)
        push!(my_node_list, curr_node)
        temp_name += 1
        Random.shuffle!(my_node_list)
    end # while
    root = pop!(my_node_list)
    set_binary!(root)

    return root
end # function create_tree_from_leaves


"""
    post_order(root::Node, traversal::Vector{Node})::Vector{Node}

This function performs a post order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the post order traversal.
"""
function post_order(root::Node, traversal::Vector{Node})::Vector{Node}
   if root.nchild != 0
       for index in 1:root.nchild
            post_order(root.child[index], traversal)
       end # for
   end # if
   push!(traversal, root)
   return traversal
end # function post_order


"""
    post_order(root:Node)::Vector{Node}

This function does post order traversal. It is meant as a wrapper. Only the root
node needs to be supplied.
"""
function post_order(root::Node)::Vector{Node}
    t::Vector{Node} = []
    post_order(root, t)
    return t
end # function post_order

"""
    tree_length(root::Node)::Float64

This function calculates the tree_length.
"""
function tree_length(root::Node)::Float64
    l::Float64 = 0.0
    for node in post_order(root)
        l += node.inc_length
    end # for
    return l
end # function tree_length


"""
    tree_height(root::Node)::Float64

This function calculates the tree height.
"""
function tree_height(root::Node)::Float64
    max_len = -Inf
    for node in post_order(root)
        if node.nchild == 0
            temp = path_length(root, node)
            if temp > max_len
                max_len = temp
            end # if
        end # end if
    end # end for
    return max_len
end # function tree_height


"""
    path_length(ancestor::Node, descendant::Node)::Float64

Note: The function assumes there is an ancestral relationship between the two nodes.

This function calculates the length of the path separating the ancestor from the
offspring node. The function follows the path specified through the binary
description of the node.
"""
function path_length(ancestor::Node, descendant::Node)::Float64
    l::Float64 = 0

    for  i in descendant.binary[length(ancestor.binary)+1:end]
        println(i)
        println(typeof(i))
        if i == '0'
            ancestor = ancestor.child[2]
        else
            ancestor = ancestor.child[1]
        end # if
        l += ancestor.inc_length
    end # for
    return l
end #function path_length


"""
    set_binary!(root::Node)

Assign a binary representation to each node, which specifies the path from the
root to this node via the binary representation of the node.
A left turn is a 1 in binary and a right turn a 0.
"""
function set_binary!(root::Node)
    if root.root == true
        root.binary = "1"
    end # if
    if root.nchild != 0
        left::Node = root.child[1]
        right::Node = root.child[2]
        left.binary = string(root.binary, "1")
        right.binary = string(root.binary, "0")
        set_binary!(left)
        set_binary!(right)
    end # if
end # function set_binary


"""
    get_leaves(root::Node)::Vector{Node}

Get all the leaves of this Node. It is meant as a wrapper, only the root node
needs to be supplied
"""
function get_leaves(root::Node)::Vector{Node}
    leave_list::Vector{Node} = []
    get_leaves(root, leave_list)
    return leave_list

end # function get_leaves


"""
    get_leaves(root::Node, leave_list::Vector{Node})::Vector{Node}

Get all the leaves of a node. All leaves of this node are visited using a
recursive application of the
`get_leaves(root::Node, leave_list::Vector{Node})::Vector{Node}` function
"""
function get_leaves(root::Node, leave_list::Vector{Node})::Vector{Node}
    if root.nchild == 0
        push!(leave_list, root)
    else
        for child in root.child
            get_leaves(child, leave_list)
        end # for
    end # if
    return leave_list
end # function get_leaves


"""
    random_node(root::Node)::Node

This function returns a random node from the tree.
"""
function random_node(root::Node)::Node
    post_order_trav = post_order(root)
    return rand(post_order_trav)
end # function random_node


"""
    move!(node1::Node, node2::Node, proportion::Float64)

Change the incomming length of node1 and node2 while keeping there combined length
constant.
"""
function move!(node1::Node, node2::Node, proportion::Float64)
    total::Float64 = node1.inc_length + node2.inc_length
    fp::Float64 = total*proportion
    sp::Float64 = total-fp

    node1.inc_length = fp
    node2.inc_length = sp

end # function move!


"""
    get_branchlength_vector(post_order::Vector{Node})::Vector{Float64}

Return a vector of branch lenghts.
"""
function get_branchlength_vector!(post_order::Vector{Node})::Vector{Float64}
    return [n.inc_length for n in post_order]
end # function get_branchlength_vector


"""
    get_branchlength_vector!(root::Node)::Vector{Float64}

Return a vector of branch lenghts.
"""
function get_branchlength_vector!(root::Node)::Vector{Float64}
    return get_branchlength_vector!(post_order(root))
end # function get_branchlength_vector


"""
    set_branchlength_vector!(post_order::Vector{Node}, blenvec::Vector{Float64})::Vector{Node}

This function sets the branch lengths of a tree to the values specified in blenvec.
"""
function set_branchlength_vector!(post_order::Vector{Node}, blenvec)::Vector{Node}
    @assert size(post_order) == size(blenvec)
    for (ind,node) in enumerate(post_order)
        node.inc_length = blenvec[ind]
    end # for
    return post_order
end # function set_branchlength_vector!


"""
    set_branchlength_vector!(root::Node, blenvec::Vector{Float64})::Node

This function sets the branch lengths of a tree to the values specified in blenvec.
"""
function set_branchlength_vector!(root::Node, blenvec::Vector{Float64})::Node
    return last(set_branchlength_vector!(post_order(root), blenvec))
end # function set_branchlength_vector!


"""
    get_sum_seperate_length!(root::Node)::Vector{Float64}

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(root::Node)::Vector{Float64}
    return get_sum_seperate_length!(post_order(root))
end # function get_sum_seperate_length!


"""
    get_sum_seperate_length!(root::Node)::Vector{Float64}

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(post_order::Vector{Node})::Vector{Float64}
    res_int::Float64 = 0.0
    res_leave::Float64 = 0.0
    for node in post_order
        if node.nchild != 0
            # internal branches
            res_int += node.inc_length
        else
            # branches leading to leaves
            res_leave += node.inc_length
        end # if
    end # for
    return [res_int, res_leave]
end # function get_sum_seperate_length!


#"""
# proper Markdown comments are not possible
#
#This part creates functions which enable the search for different nodes in the
#tree. It is possible to look for a node via its name, its binary representation
#or to find the root.
#This functionality can be extended by adding more fields to the nodes and the
#meta programmming part here.
#"""
for (sym, my_type) in [(:binary, :String), (:name, :String), (:root ,:Bool)]
    # extend the list to look for more fields in the node
    @eval function $(Symbol(string("find_by_$sym")))(tree::Node, identifier::$my_type)::Node
        # create each function and make it so it only accepts the correct type
        local all_nodes = post_order(tree) # make sure all_nodes only belongs to this function
        for node in all_nodes
            if node.$sym == identifier
                # return the node if it is found
                return node
            end # if
        end # for
        # the node is not found. Therefore throw an error!
        throw("The node identified by $identifier is not in the tree.")
    end # function
end







#end # module my_tree

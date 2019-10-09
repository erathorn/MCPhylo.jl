#=
my_tree:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#

#TODO: Automate export of automatically genereated funtions

function Base.summary(io::IO, d::Node)
    summary(io, d.name)
end

function Base.show(io::IO, d::Node)
    show(io, d.name)
end

function showall(io::IO, d::Node)
  show(io, d)
  print(io, "\nNode:\n")
  show(io, "text/plain", d.name)
  print(io, "\n\n#children:\n")
  show(io, d.nchild)
  print(io, "\n\nbinary:\n")
  show(io, d.binary)
end




"""
    add_child(mother_node::Node, child::Node)

This function adds a child to the mother node.
The arity of the mother node is increased by `1` and the root
status of the child is set to `False`.
"""
function add_child!(mother_node::Node, child::Node, left::Bool)
    if left
        mother_node.lchild = child
        child.mother = mother_node
    else
        mother_node.rchild = child
        child.mother = mother_node
    end
    #push!(mother_node.child, child)
    mother_node.nchild += 1
    child.root = false
end # function add_child

"""
    add_child(mother_node::Node, child::Node)

This function adds a child to the mother node.
The arity of the mother node is increased by `1` and the root
status of the child is set to `False`.
"""
function add_child!(mother_node::Node, child::Node)
    if ismissing(mother_node.lchild)
        mother_node.lchild = child
        child.mother = mother_node
    elseif ismissing(mother_node.rchild)
        mother_node.rchild = child
        child.mother = mother_node

    else
        throw("The node $mother_node already has two children.")
    end
    #add_child!(mother_node, child, id)
    mother_node.nchild += 1
    child.root = false
end # function add_child


"""
    remove_child!(mother_node::Node, index::int)Node

This function removes a child from the list of nodes which are daughters of this
node. The removed node is returned.
"""
function remove_child!(mother_node::Node, left::Bool)::Node
    if left
        rv = mother_node.lchild
        mother_node.lchild = missing
        rv.mother = missing
    else
        rv = mother_node.rchild
        mother_node.rchild = missing
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
    if mother_node.lchild == child
        mother_node.lchild = missing
        child.mother = missing

    elseif mother_node.rchild == child
        mother_node.rchild = missing
        child.mother = missing
    else
        # the node is not found. Therefore throw an error!
        throw("The node $child is not a child of $mother_node.")
    end # end if
    mother_node.nchild -= 1
    return child
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
        nn =  Node(node_name, zeros(Float64, (2, node_size)),missing, missing, missing, 0, true, 0.0, "0", 0)
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
        curr_node::Node = Node(string(temp_name), zeros(Float64, (2, node_size)), missing, missing, missing, 0, true, 0.0, "0", 0)
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


"""
    post_order(root::Node, traversal::Vector{Node})::Vector{Node}

This function performs a post order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the post order traversal.
"""
function post_order(root::Node, traversal::Vector{Node})::Vector{Node}
   if root.nchild != 0
        isdefined(root, :lchild) && post_order(root.lchild, traversal)
        isdefined(root, :rchild) && post_order(root.rchild, traversal)
   end # if
   push!(traversal, root)
   return traversal
end # function post_order_trav


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
    pre_order(root::Node, traversal::Vector{Node})::Vector{Node}

This function performs a pre order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the pre order traversal.
"""
function pre_order(root::Node, traversal::Vector{Node})::Vector{Node}
    push!(traversal, root)
    if root.nchild != 0
        pre_order(root.lchild, traversal)
        pre_order(root.rchild, traversal)
    end # if
    return traversal
end # function pre_order!


"""
    pre_order(root:Node)::Vector{Node}

This function does pre order traversal. It is meant as a wrapper. Only the root
node needs to be supplied.
"""
function pre_order(root::Node)::Vector{Node}
    t::Vector{Node} = []
    pre_order(root, t)
    return t
end # function pre_order!

function newick(root::Node)
    newickstring = newick(root, "")
    newickstring = string(newickstring, ";")
    return newickstring
end
function newick(root::Node, newickstring::AbstractString)
    if root.nchild != 0
        newickstring = string(newickstring, "(")
        newickstring = newick(root.lchild, newickstring)
        newickstring = string(newickstring, ",")
        newickstring = newick(root.rchild, newickstring)
        newickstring = string(newickstring, ")")
        newickstring = string(newickstring, root.name)
        newickstring = string(newickstring, ":")
        newickstring = string(newickstring, root.inc_length)
        return newickstring
    else
        newickstring = string(newickstring, root.name)
        newickstring = string(newickstring, ":")
        return string(newickstring, root.inc_length)
    end
end

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
    mn = descendant
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
function get_sister(node::Node)::Node
    mother = node.mother
    if node == mother.lchild
        return mother.rchild
    else
        return mother.lchild
    end # if
end # function


"""
    get_mother(root::Node, node::Node)::Node

This function gets the mother of `node`. It does so by looking for the respective
binary representation of the mother node.
"""
function get_mother(node::Node)::Node
    return node.mother
end # function

"""
    set_binary!(root::Node)

Assign a binary representation to each node, which specifies the path from the
root to this node via the binary representation of the node.
A left turn is a 1 in binary and a right turn a 0.
"""
function set_binary!(root::Node)
    if root.root
        root.binary = "1"
    end # if
    if root.nchild != 0
        left::Node = root.lchild
        right::Node = root.rchild
        left.binary = string(root.binary, "1")
        right.binary = string(root.binary, "0")
        set_binary!(left)
        set_binary!(right)
    end # if
end # function set_binary

"""
    number_nodes!(root::Node)::Nothing

This function assigns a unique, sequential number to each node.
"""
function number_nodes!(root::Node)::Nothing
    for (index, value) in enumerate(post_order(root))
        value.num = index
    end # for
end # fuction number_nodes


"""
    get_leaves(root::Node)::Vector{Node}

Get all the leaves of this Node. It is meant as a wrapper, only the root node
needs to be supplied
"""
function get_leaves(root::Node)::Vector{Node}
    leave_list::Vector{Node} = [i for i in post_order(root) if i.nchild == 0]
    #get_leaves(root, leave_list)
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
function get_branchlength_vector(post_order::Vector{Node})::Vector{Float64}
    out = zeros(length(post_order)-1)
    for i in eachindex(post_order)
        if !post_order[i].root
            out[post_order[i].num]= post_order[i].inc_length
        end
    end
    return out
end # function get_branchlength_vector


"""
    get_branchlength_vector!(root::Node)::Vector{Float64}

Return a vector of branch lenghts.
"""
function get_branchlength_vector(root::Node)::Vector{Float64}
    return get_branchlength_vector(post_order(root))
end # function get_branchlength_vector

function get_branchlength_vector(t::TreeStochastic)
    get_branchlength_vector(t.value)
end # function

"""
    set_branchlength_vector!(post_order::Vector{Node}, blenvec::Vector{Float64})::Vector{Node}

This function sets the branch lengths of a tree to the values specified in blenvec.
"""
function set_branchlength_vector!(post_order::Vector{Node}, blenvec::Array{Float64})::Vector{Node}
    @assert length(post_order)-1 == length(blenvec)
    for node in post_order
        if !node.root
            node.inc_length = blenvec[node.num]
        end
    end # for
    return post_order
end # function set_branchlength_vector!


"""
    set_branchlength_vector!(root::Node, blenvec::Vector{Float64})::Node

This function sets the branch lengths of a tree to the values specified in blenvec.
"""
function set_branchlength_vector!(root::Node, blenvec::Array{Float64})::Node
    return last(set_branchlength_vector!(post_order(root), blenvec))
end # function set_branchlength_vector!


function set_branchlength_vector!(t::TreeStochastic, blenvec::Array{Float64})
    set_branchlength_vector!(t.value, blenvec::Array{Float64})
end # function

function set_branchlength_vector!(t::TreeStochastic, blenvec::ArrayStochastic)
    set_branchlength_vector!(t.value, blenvec.value)
end # function


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

function internal_external_map(t::TreeStochastic)::Vector{Int64}
    internal_external_map(t.value)
end

function internal_external_map(root::Node)::Vector{Int64}
    internal_external_map(post_order(root))
end

function internal_external_map(post_order::Vector{Node})::Vector{Int64}
    my_map = Vector{Int64}(undef, length(post_order))
    for node_ind in eachindex(post_order)
        if post_order[node_ind].nchild != 0
            my_map[node_ind] = 1
        else
            my_map[node_ind] = 0
        end
    end
    return my_map
end

#"""
# proper Markdown comments are not possible
#
#This part creates functions which enable the search for different nodes in the
#tree. It is possible to look for a node via its name, its binary representation
#or to find the root.
#This functionality can be extended by adding more fields to the nodes and the
#meta programmming part here.
#"""
for (sym, my_type) in [(:binary, :String), (:name, :String), (:root ,:Bool), (:num, :Int64)]
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

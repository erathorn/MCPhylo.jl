include("../MCPhylo.jl")
using .MCPhylo
function ladderize_tree!(root::T,ascending::Bool=true)::Nothing where T<:AbstractNode
    if root.nchild == 0
        return
    end
    ndescendants = Array{Float64,1}(undef, length(root.children))
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

using Serialization

function newick(root::T)::String  where T<:Node
    # get the newickstring
    newickstring = newick(root, "")

    # Some polishing.
    newickstring = chop(newickstring)
    newickstring = string(newickstring, ";")
    return newickstring
end

function newick(root::T, newickstring::AbstractString) where T<:AbstractNode
    if root.nchild != 0
        # internal node
        newickstring = string(newickstring, "(")
        for child in root.children
            newickstring = string(newick(child,newickstring))
        end # for
        newickstring = chop(newickstring)
        return string(newickstring,")", root.name, ":", root.inc_length,",")

    else
        # leave
        return string(newickstring, root.name, ":", root.inc_length, ",")
    end # if
end

function node_height(root::T)  where T<:Node
    for node in post_order(root)
        if node.nchild == 0
            node.height = 0
        else
            node.height = maximum([child.inc_length+child.height for child in node.children])
        end
    end
end # function node_height
function post_order(root::T, traversal::Vector{T})::Vector{T} where T<:AbstractNode
   if root.nchild != 0
        for child in root.children
            post_order(child, traversal)
        end
   end # if
   push!(traversal, root)
   return traversal
end # function post_order_trav
function post_order(root::T)::Vector{T} where T<:Node
    t::Vector{T} = []
    post_order(root, t)
    return t
end # function post_order

"""
    neighbor_joining(dm::Array{Float64,2}, Array{String,1})

This function returns a phylogenetic tree by using neighbor-joining based on a
given distance matrix and an array of leaf names. Returns a node of the
resulting tree, from which it can be traversed.
"""
function neighbor_joining(dm::Array{Float64,2}, leaf_names::Array{String,1})
    n = size(dm)[1]
    if n != size(leaf_names, 1)
        throw("Distance Matrix and leaf names array size do not match")
    end # end if
    # build array of leaves from leaf names
    leaves = Array{Node,1}(undef, size(dm, 1))
    for (ind, leaf) in enumerate(leaf_names)
        new_leaf = Node(leaf)
        leaves[ind] = new_leaf
    end # end for
    neighbor_joining_int(dm, leaves)
end

"""
    neighbor_joining(dm::Array{Float64,2})

This function returns a phylogenetic tree by using neighbor-joining based on a
given distance matrix. Creates an array of nodes to be used as leaves. Returns
a node of the resulting tree, from which it can be traversed.
"""
function neighbor_joining(dm::Array{Float64,2})
    n = size(dm)[1]
    leaves = Array{Node,1}(undef, n)
    # build array of dummy leaves
    for i in 1:n
        new_leaf = Node()
        leaves[i] = new_leaf
    end # end for
    neighbor_joining_int(dm, leaves)
end

"""
    neighbor_joining_int(dm::Array{Float64,2},leaves::Vector{Node})

Internal function that is called by both neighbor_joining methods. Contains the
actual neighbor-joining algorithm, i.e. builds a phylogenetic tree from the
given distance matrix and array of leaves. Returns a node of that tree, from
which it can be traversed.
"""
function neighbor_joining_int(dm::Array{Float64,2}, leaves::Vector{Node})

    n = size(dm)[1]
    # count for node names
    count = 1
    while n > 1
        # build Q1 matrix
        q1 = zeros(Float64, n, n)
        for i in 1:n
            for j in i+1:n
                q1[i, j] =
                    q1[j, i] =
                        (n - 2) * dm[i, j] - sum(dm[i, :]) - sum(dm[j, :])
            end # end for
        end # end for
        # locate minimum of Q1
        index = findmin(q1)[2]
        first_node = leaves[index[2]]
        second_node = leaves[index[1]]
        new_node = Node("Node_$count")
        count += 1
        first_node.inc_length =
            0.5 * dm[index] +
            (1 / (2 * (n - 2))) * (sum(dm[index[2], :]) - sum(dm[index[1], :]))
        second_node.inc_length = dm[index] - first_node.inc_length
        add_child1!(new_node, first_node)
        add_child1!(new_node, second_node)
        # update array with leaves
        deleteat!(leaves, [index[2], index[1]])
        insert!(leaves, 1, new_node)
        # initalize next distance matrix
        n -= 1
        next_dm = zeros(Float64, n, n)
        # fill first row and column of next distance matrix
        j = 1
        for i in 2:n
            while j == index[1] || j == index[2]
                j += 1
            end # end while
            next_dm[i, 1] =
                next_dm[1, i] =
                    0.5 *
                    (dm[index[1], j] + dm[index[2], j] - dm[index[1], index[2]])
            j += 1
        end # end for
        # copy values of last distance matrix to finish filling the next one
        array = [1:1:n+1;]
        deleteat!(array, [index[2], index[1]])
        next_dm[2:end, 2:end] .= dm[array, array]
        # add final node to tree, if only 2 nodes are left in the list of nodes
        dm = next_dm
        if n == 2
            final_leaf = pop!(leaves)
            final_leaf.inc_length = dm[1, 2]
            # add third child to new node created previously
            add_child1!(new_node, final_leaf)
            set_binary!(new_node)
            number_nodes!(new_node)
            return new_node
        end # end if
    end # end while
end # end function neighbor_joining
"""
    add_child1!(mother_node::Node, child::Node)

This function adds a child to the mother node.
The arity of the mother node is increased by `1` and the root
status of the child is set to `False`.
"""
function add_child1!(mother_node::Node, child::Node)
    push!(mother_node.children, child)
    child.mother = mother_node
    mother_node.nchild += 1
    child.root = false
    mother_node.initialized=true
end # function add_child1!

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
There are three distance matrices and lists of leaves, which you can read into
a julia object as shown here
"""
distance_matrix = deserialize("./MCPhylo/src/Tree/dm_1.jls")
leaves_list = deserialize("./MCPhylo/src/Tree/leaves_1.jls")


"""
After creating a tree with the neighbor joining method, you can save the tree as
a newick string and save it into a file as follows
"""
tree = neighbor_joining(distance_matrix, leaves_list)
f = open("newick_output.nwk", "w")
println(f, newick(tree))
close(f)


"""
Using the dendropy python package (https://dendropy.org/) you can compare your tree
to the gold standard one in newick1.nwk (respectively for newick2, and newick3)
Apparently they can also compute a neighbor joining tree from a distance matrix,
if you find it necessary, you can also use these methods for comparison.
"""

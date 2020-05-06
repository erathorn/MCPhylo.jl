using Serialization
include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo


"""
    neighbor_joining(dm::Array{Float64,2}, Array{String,1})

This function returns a phylogenetic tree by using neighbor-joining based on a
given distance matrix and an array of leaf names. Returns a node of the
resulting tree, from which it can be traversed.
"""
function neighbor_joining(
    dm::Array{Float64,2},
    leaf_names::Array{String,1}
)
    n = size(dm)[1]
    if n != size(leaf_names, 1)
        throw("Distance Matrix and leaf names array size do not match")
    end # end if
    # build array of leaves from leaf names
    leaves = Array{Node, 1}(undef,size(dm, 1))
    for (ind,leaf) in enumerate(leaf_names)
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
    leaves = Array{Node, 1}(undef, n)
    # build array of dummy leaves
    for i = 1:n
        new_leaf = Node("leaf_$i")
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
function neighbor_joining_int(
    dm::Array{Float64,2},
    leaves::Vector{Node},
)

    n = size(dm)[1]
    # count for node names
    count = 1
    while n > 1
        # build Q1 matrix
        q1 = zeros(Float64, n, n)
        for i = 1:n
            for j = i+1:n
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
        add_child!(new_node, first_node)
        add_child!(new_node, second_node)
        # update array with leaves
        deleteat!(leaves, [index[2], index[1]])
        insert!(leaves, 1, new_node)
        n -= 1
        # initalize next distance matrix
        next_dm = zeros(Float64, n, n)
        # fill first row and column of next distance matrix
        j = 1
        for i = 2:n
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
            add_child!(new_node, final_leaf)
            set_binary!(new_node)
            number_nodes!(new_node)
            return new_node
        end # end if
    end # end while
end # end function neighbor_joining

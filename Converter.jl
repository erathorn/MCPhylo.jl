module Converter

include("./Tree_Basics.jl")
using DataFrames
using Markdown
using ..Tree_Basics: Node, post_order, set_binary!, add_child!

export to_df, from_df, to_newick

#TODO: Nexus support
#TODO: Newick parser

"""
    to_df(root::Node)::DataFrame

This function returns a matrix representation of the tree structure. The matirx
is returned as a DataFrame so that the names of the columns are the names of the
tips in the tree. The entry `df[i,j]` is the length of the edge connecting node
`i` with node `j`.
"""
function to_df(root::Node)::DataFrame

    post_order_iteration = post_order(root)

    name_list = [i.name for i in post_order_iteration]

    temp_ar = zeros(Float64, (length(post_order_iteration), length(post_order_iteration)))

    for i in post_order_iteration
        if i.nchild != 0
            ind = indexin(i.name, name_list)
            for j in i.child
                ind2 = indexin(j.name, name_list)
                temp_ar[ind[1], ind2[1]] = j.inc_length
            end # end for
        end # end if
    end # end for
    df = convert(DataFrame, temp_ar)
    names!(df, [Symbol(i) for i in name_list])
    return df
end # end function to_df


"""
    from_df(df::DataFrame)::Node

This function takes a DataFrame and turns it into a tree. It assumes a rooted
binary tree is stored in the matrix. No checks are performed.
"""
function from_df(df::DataFrame)::Node

    node_list::Vector{Node} = [Node(parse(Float64, String(i)), [0.0], Node[], 0, true, 0.0, "0") for i in names(df)]

    for (col_index, col) in enumerate(eachcol(df))
        for (row_index, entry) in enumerate(col)
            if entry != 0
                # there is a branch connecting these nodes
                node_list[col_index].inc_length = entry
                add_child!(node_list[row_index], node_list[col_index])
            end # end if
        end # end for
    end # end for
    i::Int = 0
    # find the root
    for (ind, n) in enumerate(node_list)
        if n.root == true
            i = ind
            break
        end # end if
    end # for
    node::Node = node_list[i]

    # do some bookeeping here and set the binary representation of the nodes
    set_binary!(node)
    return node
end # function from_df

function to_newick(node::Node)::String
    ret_str = ""
    if node.nchild == 0
        return string(node.name)*":"*string(node.inc_length)
    else
        ret_str *= "(" *to_newick(node.child[1])* "," *to_newick(node.child[2])*")"*string(node.name)*":"*string(node.inc_length)
    end # if
    if node.root == true
        return ret_str*";"
    else
        return ret_str
    end # if
end # function



end # module Converter

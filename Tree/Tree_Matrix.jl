



"""
    create_tree_from_leaves(leaf_nodes::Vector{String})::Array{Float64}

This function creates a  random binary tree from a list of leaf nodes.
The matrix representing the tree structure gets returned
"""
function create_tree_from_leaves_mat(leaf_nodes::Vector{String})::Array{Float64}

    l::Int64 = length(leaf_nodes)
    my_node_mat::Array{Float64, 2} = zeros((l+l-1,l+l-1))

    my_node_list = collect(1:l)
    # Internal nodes are created using integers as names.
    temp_name::Int64 = l+1

    # shuffle the node list to get a random tree
    Random.shuffle!(my_node_list)

    while length(my_node_list) != 1
        # get two nodes
        # create a new mother node to which the two first nodes are added as children
        # add the new mother node to the list and reshuffle
        first_child::Int64 = pop!(my_node_list)
        second_child::Int64 = pop!(my_node_list)

        my_node_mat[temp_name, first_child] = rand(Uniform(0,1))
        my_node_mat[temp_name, second_child] = rand(Uniform(0,1))

        push!(my_node_list, temp_name)

        temp_name += 1
        Random.shuffle!(my_node_list)
    end # while

    return my_node_mat
end # function create_tree_from_leaves

"""
    function get_neighbours(vec::Array{Float64,1})::Array{Int64}

This function returns an array containing the neighbours of the node. If a row
is supplied as argument, the daughter nodes are returned, otherwise the mother node.

The function is equiped the with @inline decorator, since it will be used quite
frequently. This should improve speed.
"""
@inline function get_neighbours(vec::Array{Float64,1})::Array{Int64}
    l::Int64 = size(vec)[1]
    r::Vector{Int64} = []
    @inbounds for i in 1:l
        if vec[i] > 0.0
            push!(r, i)
        end # if
    end # for
    return r
end # get_neighbours

"""
    find_root(mat::Array{Float64,2})::Int64

This function returns the root node of the tree speciefied through the matrix.
"""
function find_root(mat::Array{Float64,2})::Int64
    l::Int64 = size(mat)[1]
    for i in 1:l
        if length(get_neighbours(mat[:,i])) == 0
            return i
        end # if
    end # for
    return ErrorException("No root node found")
end # find_root


"""
    post_order(mat::Array{Float64,2}, node::Int64, traversal::Array{Int64, 1})::Array{Int64}

This function performs a post order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the post order traversal.
"""
function post_order(mat::Array{Float64,2}, node::Int64, traversal::Vector{Int64})::Vector{Int64}

    for i in get_neighbours(mat[node, :])
        post_order(mat, i, traversal)
    end # for
    push!(traversal, node)
    return traversal
end # post_order


"""
    post_order(mat::Array{Float64, 2})::Array{Int64}

This function performs a post order traversal through the tree.
"""
function post_order(mat::Array{Float64, 2})::Vector{Int64}
    root::Int64 = find_root(mat)
    traversal::Vector{Int64} = []
    return post_order(mat, root, traversal)
end # post_order


"""
    pre_order(mat::Array{Float64,2}, node::Int64, traversal::Array{Int64, 1})::Array{Int64}

This function performs a pre order traversal through the tree. It is assumed that `root` is the
root of the tree. Thus, if `root` is not the root, the subtree defined by the root `root` is
used for the pre order traversal.
"""
function pre_order(mat::Array{Float64,2}, node::Int64, traversal::Vector{Int64})::Vector{Int64}
    push!(traversal, node)
    for i in get_neighbours(mat[node, :])
        post_order(mat, i, traversal)
    end # for
    return traversal
end # pre_order


"""
    pre_order(mat::Array{Float64, 2})::Array{Int64}

This function performs a post order traversal through the tree.
"""
function pre_order(mat::Array{Float64, 2})::Vector{Int64}
    root::Int64 = find_root(mat)
    traversal::Vector{Int64} = []
    return post_order(mat, root, traversal)
end # pre_order


"""
    tree_length(mat::Array{Float64,2}) = sum(mat)

Return the entire tree length.
"""
tree_length(mat::Array{Float64,2}) = sum(mat)


"""
    get_leaves(mat::Array{Float64,2})::Array{Int64,1}

get all the leaves of the tree specified by mat.
"""
function get_leaves(mat::Array{Float64,2})::Vector{Int64}
    l::Int64 = size(mat)[1]
    leaves::Vector{Int64} = []
    for i in 1:l
        if length(get_neighbours(mat[i,:])) == 0
            push!(leaves, i)
        end # if
    end # for
    return leaves
end # get_leaves


"""
    get_mother(mat::Array{Float64,2}, node::Int64)::Int64

get the mother node of the node supplied as an argument in the tree specied by mat.
"""
function get_mother(mat::Array{Float64,2}, node::Int64)::Int64
    l::Int64 = size(mat)[1]
    for i in 1:l
        if mat[i,node] > 0
            return i
        end # if
    end # for
    return ErrorException("No mother node found for node $node")
end # get_mother


"""
    get_sum_seperate_length!(mat::Array{Float64,2})::Vector{Float64}

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(mat::Array{Float64,2})::Vector{Float64}
    leaves::Vector{Int64} = get_leaves(mat)
    l::Int64 = size(mat)[1]
    res_int::Float64 = 0.0
    res_leave::Float64 = 0.0
    for i in 1:l
        if i in leaves
            # branches leading to leaves
            res_leave += sum(mat[:,i])
        else
            # internal branches
            res_int += sum(mat[:,i])
        end # if
    end # for
    return [res_int, res_leave]
end # function get_sum_seperate_length!


"""
    make_tree_with_data_mat(filename::String)

This function creates a tree where the terminal nodes get the data specified in
the NEXUS file.
"""
function make_tree_with_data_mat(filename::String)
    # get all the information from the NEXUS file
    n_tax, nc, gap, miss, df = ParseNexus(filename)

    # create random tree
    new_tree = create_tree_from_leaves_mat(df[:Language])
    l = size(new_tree)[1]
    data_arr::Array{Float64,3} = zeros(2, nc, l)

    for (lind, row) in enumerate(eachrow(df))

        for (ind, i) in enumerate(row.Data)
            if i == '0'
                data_arr[1,ind,lind] = 1.0
            elseif i == '1'
                data_arr[2,ind, lind] = 1.0
            else
                data_arr[1, ind, lind] = 1.0
                data_arr[2, ind, lind] = 1.0
            end # if
        end # for
    end
    return new_tree, data_arr, df
end # function make_tree_with_data

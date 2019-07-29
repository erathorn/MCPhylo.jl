



"""
    create_tree_from_leaves(leaf_nodes::Vector{T})::Node

This function creates a  random binary tree from a list of leaf nodes.
The root node as access point for the tree is returned.
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


@inline function get_neighbours(vec::Array{Float64,1})::Array{Int64}
    l::Int64 = size(vec)[1]
    r::Array{Float64} = []
    @inbounds for i in 1:l
        if vec[i] > 0.0
            push!(r, i)
        end
    end
    return r
end

function find_root(mat::Array{Float64,2})::Int64
    l::Int64 = size(mat)[1]
    for i in 1:l
        if length(get_neighbours(mat[:,i])) == 0
            return i
        end
    end
    return ErrorException("No root node found")
end


function post_order(mat::Array{Float64,2}, node::Int64, traversal::Array{Int64, 1})::Array{Int64}

    for i in get_neighbours(mat[node, :])
        post_order(mat, i, traversal)
    end
    push!(traversal, node)
    return traversal
end

function post_order(mat::Array{Float64, 2})::Array{Int64}
    root::Int64 = find_root(mat)
    traversal::Array{Int64,1} = []
    return post_order(mat, root, traversal)
end

function pre_order(mat::Array{Float64,2}, node::Int64, traversal::Array{Int64, 1})::Array{Int64}
    push!(traversal, node)
    for i in get_neighbours(mat[node, :])
        post_order(mat, i, traversal)
    end
    return traversal
end

function pre_order(mat::Array{Float64, 2})::Array{Int64}
    root::Int64 = find_root(mat)
    traversal::Array{Int64,1} = []
    return post_order(mat, root, traversal)
end

tree_length(mat::Array{Float64,2}) = sum(mat)

function get_leaves(mat::Array{Float64,2})::Array{Int64,1}
    l::Int64 = size(mat)[1]
    leaves::Array{Int64, 1} = []
    for i in 1:l
        if length(get_neighbours(mat[i,:])) == 0
            push!(leaves, i)
        end
    end
    return leaves
end

function get_mother(mat::Array{Float64,2}, node::Int64)::Int64
    l::Int64 = size(mat)[1]
    for i in 1:l
        if mat[i,node] > 0
            return i
        end
    end
    return ErrorException("No mother node found for node $node")
end

function get_sum_seperate_length!(mat::Array{Float64,2})::Vector{Float64}
    leaves::Array{Int64} = get_leaves(mat)
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


function make_tree_with_data_mat(filename::String)#::Tree_Module.Node
    # get all the information from the NEXUS file
    n_tax, nc, gap, miss, df = ParseNexus(filename)

    # create random tree
    new_tree = create_tree_from_leaves_mat(df[:Language])
    l = size(new_tree)[1]
    data_arr::Array{Float64,3} = zeros(n_tax, 2, nc)

    for (lind, row) in enumerate(eachrow(df))

        for (ind, i) in enumerate(row.Data)
            if i == '0'
                data_arr[lind,1,ind] = 1.0
            elseif i == '1'
                data_arr[lind,2,ind] = 1.0
            else
                data_arr[lind,1, ind] = 1.0
                data_arr[lind,2, ind] = 1.0
            end # if
        end # for
    end
    return new_tree, data_arr, df
end # function make_tree_with_data

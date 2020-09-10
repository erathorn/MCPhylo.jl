"""
    find_common_clusters(ref_tree, tree:T)
        ::Dict{Node, Bool} where T<:AbstractNode

Use Day's algorithm to create a dictionary, that tells us for each node of the
second input tree, if its corresponding cluster is a common cluster of the two trees
"""
function find_common_clusters(ref_tree::T, tree::T)::Dict{Node, Bool} where T<:AbstractNode
    ref_leaf_names = [leaf.name for leaf in get_leaves(ref_tree)]
    leaf_names = [leaf.name for leaf in get_leaves(tree)]
    if length(ref_leaf_names) != length(leaf_names)
        throw(ArgumentError("The leaf label sets of the trees are not identical"))
    elseif length(ref_leaf_names) != length(Set(append!(ref_leaf_names, leaf_names)))
        throw(ArgumentError("The leaf label sets of the trees are not identical"))
    end
    nodes = post_order(ref_tree)
    leaves_dict = Dict{String, Int64}()
    clusters = Dict{Node, Vector{String}}()
    cluster_dict = Dict{Tuple{Int64, Int64}, Int64}()
    dict_count = 1
    for node in nodes
        if node.nchild == 0
            leaves_dict[node.name] = dict_count
            dict_count += 1
        else
            cluster = Vector{String}()
            for child in node.children
                child.nchild == 0 ? push!(cluster, child.name) : append!(cluster, clusters[child])
            end
            clusters[node] = cluster
            if node.root == true || node == node.mother.children[1]
                last_index = leaves_dict[last(cluster)]
                first_index = leaves_dict[first(cluster)]
                cluster_dict[(first_index, last_index)] = last_index
            else
                last_index = leaves_dict[last(cluster)]
                first_index = leaves_dict[first(cluster)]
                cluster_dict[(first_index, last_index)] = first_index
            end # if
        end # if
    end # for
    is_common_cluster = Dict{Node, Bool}()
    clusters = Dict{Node, Vector{String}}()
    nodes = post_order(tree)
    for node in nodes
        if node.nchild == 0
            is_common_cluster[node] = true
        else
            cluster = Vector{String}()
            for child in node.children
                child.nchild == 0 ? push!(cluster, child.name) : append!(cluster, clusters[child])
            end
            clusters[node] = cluster
            cluster_indeces = [leaves_dict[leaf] for leaf in cluster]
            if length(cluster) != maximum(cluster_indeces) - minimum(cluster_indeces) + 1
                is_common_cluster[node] = false
            elseif haskey(cluster_dict, (minimum(cluster_indeces), maximum(cluster_indeces)))
                is_common_cluster[node] = true
            else
                is_common_cluster[node] = false
            end # if/else
        end # if/else
    end # for
    return is_common_cluster
end


"""
    one_way_compatible(ref_tree::T, tree::T)::T where T<:AbstractNode

Takes two trees and returns a copy of the first one, where all the clusters that
are not compatible with the second tree are removed.
"""
function one_way_compatible(ref_tree::T, tree::T)::T where T<:AbstractNode
    ref_tree_copy = deepcopy(ref_tree)
    ref_nodes = post_order(ref_tree_copy)
    node_leafs = Dict{Node, Int64}()
    leaves_dict = Dict{String, Int64}()
    dict_count = 1
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = dict_count
            dict_count += 1
        else
            leaves_count = 0
            for child in ref_node.children
                child.nchild == 0 ? leaves_count += 1 : leaves_count += node_leafs[child]
            end # for
            node_leafs[ref_node] = leaves_count
        end # if/else
    end # for
    nodes = post_order(tree)
    leaves = Vector{Node}()
    cluster_start_indeces = Dict{Node, Int64}()
    node_binaries = Dict{String, Node}()
    for node in nodes
        if node.nchild != 0
            cluster_start_indeces[node] = cluster_start_indeces[node.children[1]]
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
        end # if / else
        node_binaries[node.binary] = node
    end # for
    leaves = order_tree!(tree, cluster_start_indeces)
    leaf_ranks = Dict(enumerate([leaf.name for leaf in leaves]))
    leaf_ranks_reverse = Dict(value => key for (key, value) in leaf_ranks)
    marked_nodes = Dict{Int64, Bool}()
    for ref_node in ref_nodes
        if ref_node.nchild != 0
            start = min_leaf_rank(leaf_ranks_reverse, ref_node)
            stop = max_leaf_rank(leaf_ranks_reverse, ref_node)
            start_node = leaves[start]
            stop_node = leaves[stop]
            lca_start_stop = node_binaries[lcp(start_node.binary, stop_node.binary)]
            length(x_left(start_node).binary) > length(lca_start_stop.binary) ?
                d = start_node : d = lca_start_stop.children[1]
            length(x_right(stop_node).binary) > length(lca_start_stop.binary) ?
                e = stop_node : e = lca_start_stop.children[end]
            if d.mother == lca_start_stop && e.mother == lca_start_stop && node_leafs[ref_node] == stop - start + 1
                marked_nodes[ref_node.num] = false
            else
                marked_nodes[ref_node.num] = true
            end # if / else
        end # if
    end # for
    for node in level_order(ref_tree_copy)
        if node.nchild != 0 && marked_nodes[node.num] == true
            delete_node!(node)
        end # if
    end # for
    return ref_tree_copy
end


"""
    merge_trees!(ref_tree::T, tree::T)::Tuple{T, Vector{T}} where T<:AbstractNode

Merge two compatible trees
"""
function merge_trees!(ref_tree::T, tree::T)::Tuple{T, Vector{T}} where T<:AbstractNode
    ref_nodes = post_order(ref_tree)
    inserted_nodes = Vector{Node}()
    leaves_dict = Dict{String, Int64}()
    dict_count = 1
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = dict_count
            dict_count += 1
        end # if
    end # for
    nodes = post_order(tree)
    ref_leaves = Vector{Node}()
    cluster_start_indeces = Dict{Node, Int64}()
    node_binaries = Dict{String, Node}()
    for node in nodes
        if node.nchild != 0
            for leaf in ref_leaves
                if node.binary == leaf.binary[1:length(node.binary)]
                    cluster_start_indeces[node] = leaves_dict[leaf.name]
                    break
                end # if
            end # for
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
            push!(ref_leaves, node)
        end # if / else
        node_binaries[node.binary] = node
    end # for
    leaves = order_tree!(tree, cluster_start_indeces)
    leaf_ranks = Dict(enumerate([leaf.name for leaf in leaves]))
    leaf_ranks_reverse = Dict(value => key for (key, value) in leaf_ranks)
    for ref_node in ref_nodes
        if ref_node.nchild != 0
            start = min_leaf_rank(leaf_ranks_reverse, ref_node)
            stop = max_leaf_rank(leaf_ranks_reverse, ref_node)
            start_node = leaves[start]
            stop_node = leaves[stop]
            lca_start_stop = node_binaries[lcp(start_node.binary, stop_node.binary)]
            xleft = x_left(start_node)
            xright = x_right(stop_node)
            length(xleft.binary) > length(lca_start_stop.binary) ?
                d = xleft : d = lca_start_stop.children[1]
            length(xright.binary) > length(lca_start_stop.binary) ?
                e = xright : e = lca_start_stop.children[end]
            if !(d == lca_start_stop.children[1] && e == lca_start_stop.children[end])
                index_d = findfirst(x -> x == d, lca_start_stop.children)
                index_e = findfirst(x -> x == e, lca_start_stop.children)
                inserted_node =
                    insert_node!(lca_start_stop, lca_start_stop.children[index_d:index_e])
                push!(inserted_nodes, inserted_node)
                set_binary!(tree)
                number_nodes!(tree)
            end # if
        end # if
    end # for
    return tree, inserted_nodes
end # function merge_trees!


"""
    order_tree!(root::T, cluster_start_indeces::Dict{T, Int64}, leaves=Vector{String}())
        ::Vector{String} where T<:AbstractNode

Helper function to order a tree based on cluster indeces
"""
function order_tree!(root::T, cluster_start_indeces::Dict{T, Int64}, leaves=Vector{T}())::Vector{T} where T<:AbstractNode
    sort!(root.children, by = child -> cluster_start_indeces[child])
    for child in root.children
        if child.nchild == 0
            push!(leaves, child)
        else
            order_tree!(child, cluster_start_indeces, leaves)
        end # if/else
    end # for
    return leaves
end


"""
    min_leaf_rank(leaf_ranks::Dict{String, Int64}, node::T)
        ::Int64 where T <: AbstractNode

Recursive helper function to find the lowest ranked leaf descendant of a node
"""
function min_leaf_rank(leaf_ranks::Dict{String, Int64}, node::T)::Int64 where T <: AbstractNode
    if node.nchild == 0
        return leaf_ranks[node.name]
    else
        possible_minima = Vector{Int64}()
        for child in node.children
            if child.nchild == 0
                push!(possible_minima, leaf_ranks[child.name])
            else
                push!(possible_minima, min_leaf_rank(leaf_ranks, child))
            end # if/else
        end # for
    end # if/else
    minimum(possible_minima)
end


"""
    max_leaf_rank(leaf_ranks::Dict{String, Int64}, node::T)
        ::Int64 where T<:AbstractNode

Recursive helper function to find the highest ranked leaf descendant of a node
"""
function max_leaf_rank(leaf_ranks::Dict{String, Int64}, node::T)::Int64 where T<:AbstractNode
    if node.nchild == 0
        return leaf_ranks[node.name]
    else
        possible_maxima = Vector{Int64}()
        for child in node.children
            if child.nchild == 0
                push!(possible_maxima, leaf_ranks[child.name])
            else
                push!(possible_maxima, max_leaf_rank(leaf_ranks, child))
            end # if/else
        end # for
    end # if/else
    maximum(possible_maxima)
end


"""
    x_left(node::T)::Tuple{T,T} where T<:AbstractNode

Helper function to find ancestor of a leaf that has said leaf as leftmost descendant
"""
function x_left(node::T)::T where T<:AbstractNode
    while true
        if node.root
            return node
        else
            mother = node.mother
            if mother.children[1] != node
                return node
            end # if
        node = mother
        end # if/else
    end # while
end


"""
    x_right(node::T)::Tuple{T,T} where T<:AbstractNode

Helper function to find ancestor of a leaf that has said leaf as rightmost descendant
"""
function x_right(node::T)::T where T<:AbstractNode
    while true
        if node.root
            return node
        else
            mother = node.mother
            if mother.children[end] != node
                return node
            end # if
        node = mother
        end # if/else
    end # while
end


"""
    majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode

Construct the majority rule consensus tree from a set of trees
"""
function majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode
    first_tree = deepcopy(trees[1])
    nodes = post_order(first_tree)
    # save leaf ranks to order the resulting tree later
    leaf_ranks = Dict{String, Int64}()
    count = 0
    for node in nodes
        if node.nchild == 0
            leaf_ranks[node.name] = count
            count += 1
        end
    end

    node_counts = convert(Vector{Int64}, ones(length(nodes)))
    count_dict = Dict(zip(nodes, node_counts))
    for tree in trees[2:end]
        nodes = level_order(first_tree)
        # delete clusters in the first tree that are not in the second
        is_common_cluster = find_common_clusters(tree, first_tree)
        for node in nodes
            if is_common_cluster[node] == true
                count_dict[node] += 1
            else
                count_dict[node] -= 1
                if count_dict[node] == 0
                    delete_node!(node)
                end # if
            end # else
        end # for

        set_binary!(first_tree)
        number_nodes!(first_tree)
        compatible_tree = one_way_compatible(tree, first_tree)
        set_binary!(compatible_tree)
        number_nodes!(compatible_tree)
        first_tree, inserted_nodes = merge_trees!(compatible_tree, first_tree)
        set_binary!(first_tree)
        number_nodes!(first_tree)
        for node in inserted_nodes
            count_dict[node] = 1
        end # for
    end # for
    set_binary!(first_tree)
    number_nodes!(first_tree)
    nodes = level_order(first_tree)
    node_counts = convert(Vector{Int64}, zeros(length(nodes)))
    count_dict = Dict(zip(nodes, node_counts))
    for tree in trees
        is_common_cluster = find_common_clusters(tree, first_tree)
        for node in nodes
            if is_common_cluster[node] == true
                count_dict[node] += 1
            end # if
        end # for
    end # for
    # delete non-majority clusters
    half = length(trees) / 2
    for node in nodes
        if count_dict[node] <= half
            delete_node!(node)
        end # if
    end # for
    # order the resulting tree
    nodes = post_order(first_tree)
    cluster_start_indeces = Dict{Node, Int64}()
    for node in nodes
        if node.nchild != 0
            cluster_start_indeces[node] = cluster_start_indeces[node.children[1]]
        else
            cluster_start_indeces[node] = leaf_ranks[node.name]
        end # if / else
    end # for
    order_tree!(first_tree, cluster_start_indeces)
    return first_tree
end

"""
    find_common_clusters(ref_tree, tree:T)
        ::Tuple{Vector{Vector{Node}}, Vector{Vector{Node}}} where T<:AbstractNode

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
    cluster_dict = Dict{Tuple{Int64, Int64}, Int64}()
    dict_count = 1
    for node in nodes
        if node.nchild == 0
            leaves_dict[node.name] = dict_count
            dict_count += 1
        else
            leaves = get_leaves(node)
            if node.root == true || node == node.mother.children[1]
                last_index = leaves_dict[last(leaves).name]
                first_index = leaves_dict[first(leaves).name]
                cluster_dict[(first_index, last_index)] = last_index
            else
                last_index = leaves_dict[last(leaves).name]
                first_index = leaves_dict[first(leaves).name]
                cluster_dict[(first_index, last_index)] = first_index
            end # if
        end # if
    end # for
    is_common_cluster = Dict{Node, Bool}()
    nodes = post_order(tree)
    for node in nodes
        if node.nchild != 0
            leaves = get_leaves(node)
            cluster_indeces = [leaves_dict[leaf.name] for leaf in leaves]
            if length(get_leaves(node)) != maximum(cluster_indeces) -
                                           minimum(cluster_indeces) + 1
                is_common_cluster[node] = false
            elseif haskey(cluster_dict, (minimum(cluster_indeces), maximum(cluster_indeces)))
                is_common_cluster[node] = true
            else
                is_common_cluster[node] = false
            end # if
        end # if
    end # for
    return is_common_cluster
end


"""
    are_compatible(cluster_names::Vector{String}, tree::T) where T<:AbstractNode

Check if a cluster of nodes (identified by their names) is compatible with a tree
"""
function are_compatible(cluster::Vector{String}, tree::T)::Bool where T<:AbstractNode
    lca = Node()
    try
        lca = find_lca(tree, cluster)
    catch e
        throw(ArgumentError("The leafs of the input leafs need to be in the tree"))
    end
    for child in lca.children
        leaves = [leaf.name for leaf in get_leaves(child)]
        if !(length(intersect(leaves, cluster)) in (0, length(leaves)))
            return false
        end # if
    end # for
    return true
end


"""
    are_compatible(cluster::Vector{T}, tree::T)::Bool where T<:AbstractNode

Check if a cluster of nodes is compatible with a tree
"""
function are_compatible(cluster::Vector{T}, tree::T)::Bool where T<:AbstractNode
    cluster = [node.name for node in cluster]
    are_compatible(cluster, tree)
end


"""
    merge_trees(ref_tree::T, tree::T)::Tuple{T, Vector{T}} where T<:AbstractNode

Merge two compatible trees
"""
function merge_trees(ref_tree::T, tree::T)::Tuple{T, Vector{T}} where T<:AbstractNode
    ref_nodes = post_order(ref_tree)
    inserted_nodes = Vector{Node}()
    leaves_dict = Dict{String, Int64}()
    dict_count = 1
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = dict_count
            dict_count += 1
        end #if
    end #for
    nodes = post_order(tree)
    cluster_start_indeces = Dict{Node, Int64}()
    node_binaries = Dict{String, Node}()
    for node in nodes
        if node.child != 0
            leaves = get_leaves(node)
            cluster_start_indeces[node] = leaves_dict[first(leaves).name]
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
        end # if
        node_binaries[node.binary] = node
    end # for
    leaves = order_tree!(tree, cluster_start_indeces)
    leaf_ranks = Dict(enumerate([leaf.name for leaf in leaves]))
    for ref_node in ref_nodes
        if ref_node.nchild != 0
            start = min_leaf_rank(leaf_ranks, ref_node)
            stop = max_leaf_rank(leaf_ranks, ref_node)
            start_node = leaf_ranks[start]
            stop_node = leaf_ranks[stop]
            lca_start_stop = node_binaries[lcp(start_node.binary, stop_node.binary)]
            length(x_left(start_node).binary) > length(lca_start_stop.binary) ?
                d = start_node : d = lca_start_stop.children[1]
            length(x_right(stop_node).binary) > length(lca_start_stop.binary) ?
                e = stop_node : e = lca_start_stop.children[end]
        end
        if d == start_node || e == stop_node
            index_d = findfirst(lca_start_stop.children, d)
            index_e = findfirst(lca_start_stop.children, e)
            inserted_node =
                insert_node!(lca_start_stop, lca_start_stop.children[index_d, index_e])
            push!(inserted_nodes, inserted_node)
            lca_start_stop.children[end].binary =
                string(lca_start_stop.binary, length(lca_start_stop.binary) - 1)
        end
    end
    return tree, inserted_nodes
end


"""
    one_way_compatible(ref_tree::T, tree::T)::T where T<:AbstractNode

Takes two trees and returns a copy of the first one, where all the clusters that
are not compatible with the second tree are removed.
"""
function one_way_compatible(ref_tree::T, tree::T)::T where T<:AbstractNode
    ref_tree_copy = deepcopy(ref_tree)
    ref_nodes = post_order(ref_tree_copy)
    leaves_dict = Dict{String, Int64}()
    dict_count = 1
    node_binaries = Dict{String, Node}()
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = dict_count
            dict_count += 1
        end #if
        if node.child != 0
            leaves = get_leaves(node)
            node_children[node] = length(leaves)
        end # if
    end #for
    nodes = post_order(tree)
    cluster_start_indeces = Dict{Node, Int64}()
    for node in nodes
        if node.child != 0
            leaves = get_leaves(node)
            cluster_start_indeces[node] = leaves_dict[first(leaves).name]
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
        end # if
        node_binaries[node.binary] = node
    end # for
    leaves = order_tree!(tree, cluster_start_indeces)
    leaf_ranks = Dict(enumerate([leaf.name for leaf in leaves]))
    marked_nodes = Dict{Int64, Bool}()
    for ref_node in ref_nodes
        if ref_node.nchild != 0
            start = min_leaf_rank(leaf_ranks, ref_node)
            stop = max_leaf_rank(leaf_ranks, ref_node)
            start_node = leaf_ranks[start]
            stop_node = leaf_ranks[stop]
            lca_start_stop = node_binaries[lcp(start_node.binary, stop_node.binary)]
            length(x_left(start_node).binary) > length(lca_start_stop.binary) ?
                d = start_node : d = lca_start_stop.children[1]
            length(x_right(stop_node).binary) > length(lca_start_stop.binary) ?
                e = stop_node : e = lca_start_stop.children[end]
            if d == start_node || e == stop_node
                if d.mother == lca_start_stop && e.mother == lca_start_stop && node_children[node] == stop - start + 1
                    marked_nodes[node.num] = true
                else marked_nodes[node.num] = false
                end # if
            end # if
        end # if
    end # for
    return delete_marked_nodes!(ref_tree_copy, marked_nodes)
end


"""
    delete_marked_nodes!(root::T, marked_nodes::Dict{Int64, Bool})::T where T<:AbstractNode

Helper function to delete marked incompatible nodes from a tree
"""
function delete_marked_nodes!(root::T, marked_nodes::Dict{Int64, Bool}) where T<:AbstractNode
    if marked_nodes[root.num] == true
        delete_node!(root)
    end
    if root.nchild != 0
        for child in root.children
            delete_marked_nodes!(child, marked_nodes)
        end # for
    end # if
end


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
        end
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
            end # if
        end # for
    end # if
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
            end # if
        end # for
    end # if
    maximum(possible_maxima)
end


"""
    x_left(node::T)::Tuple{T,T} where T<:AbstractNode

Helper function to find ancestor of node that has said leaf as leftmost descendant
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
        end # if
    end # while
end


"""
    x_right(node::T)::Tuple{T,T} where T<:AbstractNode

Helper function to find ancestor of node that has said leaf as rightmost descendant
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
        end # if
    end # while
end


"""
    majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode

Construct the majority rule consensus tree from a set of trees
"""
function majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode
    """
    all_leaves = [get_leaves(root) for root in trees]
    leaf_names = sort([leaf.name for leaf in all_leaves[1]])
    for leaves in all_leaves[2:end]
        if leaf_names != sort([leaf.name for leaf in leaves])
            throw(ArgumentError("Input Trees need to have identical leaf names"))
        end
    end
    """
    ref_tree = trees[1]
    nodes = level_order(ref_tree)
    node_counts = convert(Vector{Int64}, ones(length(nodes)))
    count_dict = Dict(zip(nodes, node_counts))
    for tree in trees[2:end]
        is_common_cluster  = find_common_clusters(tree, ref_tree)
        for node in nodes
            if is_common_clusters[node] == true
                count_dict[node] += 1
            elseif node.nchild != 0
                count_dict[node] -= 1
                if count_dict[node] == 0
                    delete_node!(node)
                end # if
            end # else
        end # for
        compatible_tree = one_way_compatible(tree, ref_tree)
        ref_tree, inserted_nodes = merge_trees!(compatible_tree, ref_tree)
        for node in inserted_nodes
            count_dict[node] = 1
        end # for
    end # for
    replace!(kv -> kv[1] => 0, count_dict)
    nodes = keys(count_dict)
    for tree in trees
        is_common_cluster = find_common_clusters(tree, ref_tree)
        if is_common_cluster[node] == true
            count_dict[node] += 1
        end # if
    end # for
    marked_nodes = Dict{Int64, Bool}
    half = length(trees) / 2
    for node in nodes
        if count_dict[node] <= half
            marked_nodes[node] = true
        else
            marked_nodes[node] = false
        end
    delete_marked_nodes!(ref_tree, marked_nodes)
    end
    return ref_tree
end

"""
    find_common_clusters(ref_tree, tree:T)
        ::Dict{Int64, Tuple{Bool, Union{Float64, Missing}}}

Use Day's algorithm to create a dictionary, that tells us for each node of the
second input tree, if its corresponding cluster is a common cluster of the trees
"""
function find_common_clusters(ref_tree::T, tree::T)::Dict{Int64, Tuple{Bool, Union{Float64, Missing}}} where T<:AbstractNode
    ref_nodes = post_order(ref_tree)
    leaves_dict = Dict{String, Tuple{Int64, Float64}}()
    clusters = Dict{Node, Vector{String}}()
    cluster_dict = Dict{Tuple{Int64, Int64}, Int64}()
    dict_count = 1
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = (dict_count, ref_node.inc_length)
            dict_count += 1
        else
            cluster = Vector{String}()
            for child in ref_node.children
                child.nchild == 0 ? push!(cluster, child.name) : append!(cluster, clusters[child])
            end # for
            clusters[ref_node] = cluster
            if ref_node.root == true || ref_node == ref_node.mother.children[1]
                last_index = leaves_dict[last(cluster)][1]
                first_index = leaves_dict[first(cluster)][1]
                cluster_dict[(first_index, last_index)] = last_index
            else
                last_index = leaves_dict[last(cluster)][1]
                first_index = leaves_dict[first(cluster)][1]
                cluster_dict[(first_index, last_index)] = first_index
            end # if/else
        end # if/else
    end # for
    for value in values(clusters)
        sort!(value)
    end # for
    clusters_reverse = Dict(value => key.inc_length for (key, value) in clusters)
    is_common_cluster = Dict{Int64, Tuple{Bool, Union{Float64, Missing}}}()
    clusters = Dict{Int64, Vector{String}}()
    leaf_count = 0
    nodes = post_order(tree)
    for node in nodes
        if node.nchild == 0
            leaf_count += 1
            try
                leaves_dict[node.name]
            catch KeyError
                throw(ArgumentError("The leafs sets of the trees need to be identical"))
            end # try
            is_common_cluster[node.num] = (true, leaves_dict[node.name][2])
        else
            cluster = Vector{String}()
            for child in node.children
                child.nchild == 0 ? push!(cluster, child.name) : append!(cluster, clusters[child.num])
            end # for
            clusters[node.num] = cluster
            cluster_indeces = [leaves_dict[leaf][1] for leaf in cluster]
            if length(cluster) != maximum(cluster_indeces) - minimum(cluster_indeces) + 1
                is_common_cluster[node.num] = (false, missing)
            elseif haskey(cluster_dict, (minimum(cluster_indeces), maximum(cluster_indeces)))
                is_common_cluster[node.num] = (true, clusters_reverse[sort!(cluster)])
            else
                is_common_cluster[node.num] = (false, missing)
            end # if/else
        end # if/else
    end # for
    if leaf_count != length(keys(leaves_dict))
        throw(ArgumentError("The leafs sets of the trees need to be identical"))
    end # if
    return is_common_cluster
end # function find_common_clusters


"""
    one_way_compatible(ref_tree::T, tree::T)::T where T<:AbstractNode

Takes two trees and returns a copy of the first one, where all the clusters that
are not compatible with the second tree are removed.
"""
function one_way_compatible(ref_tree::T, tree::T)::T where T<:AbstractNode
    # First Phase -> move to own function since it is equivalent to first phase of merge_trees
    set_binary!(tree)
    ref_tree_copy = deepcopy(ref_tree)
    ref_nodes = post_order(ref_tree_copy)
    leaves_dict = Dict{String, Int64}()
    count = 1
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = count
            count += 1
        end # if/else
    end # for
    nodes = post_order(tree)

    cluster_start_indeces = Dict{Node, Int64}()
    for node in nodes
        if node.nchild != 0
            cluster_start_indeces[node] = minimum(cluster_start_indeces[n] for n in node.children)
            #cluster_start_indeces[node.children[1]]
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
        end # if / else
    end # for
    leaves::Vector{Node} = order_tree!(tree, cluster_start_indeces)

    xleft_dict, xright_dict = depth_dicts(leaves)

    leaf_ranks_reverse = Dict(node.name => ind for (ind, node) in enumerate(leaves))

    marked_nodes = Dict{Int64, Bool}()
    for ref_node in ref_nodes
        if ref_node.nchild != 0
            start = min_leaf_rank(leaf_ranks_reverse, ref_node)
            stop = max_leaf_rank(leaf_ranks_reverse, ref_node)
            start_node = leaves[start]
            stop_node = leaves[stop]
            xleft, left_path = xleft_dict[start_node]
            xright, right_path = xright_dict[stop_node]

            if intersect(values(left_path), values(right_path)) == []
                marked_nodes[ref_node.num] = true
                continue
            end
            !xleft.root && delete!(left_path, length(split(xleft.mother.binary,",")))
            !xright.root && delete!(right_path, length(split(xright.mother.binary,",")))
            depth_left = length(split(xleft.binary,","))
            depth_right = length(split(xright.binary, ","))
            depth = depth_left >= depth_right
            depth ? p2 = right_path[depth_left] : p2 = left_path[depth_right]
            depth ? p1 = xleft : p1 = xright
            if p2 == p1
                r = p1
                if depth
                    d = left_path[depth_left + 1]
                    e = right_path[depth_left + 1]
                else
                    d = left_path[depth_right + 1]
                    e = right_path[depth_right + 1]
                end # if/else
            elseif p1.mother == p2.mother
                r = p1.mother
                depth ? d = p1 : d = p2
                depth ? e = p2 : e = p1
            else
                marked_nodes[ref_node.num] = true
                continue
            end # if/else
            if length(d.mother.binary) <= length(r.binary) && length(e.mother.binary) <= length(r.binary)
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
    set_binary!(ref_tree_copy)
    number_nodes!(ref_tree_copy)
    return ref_tree_copy
end # function one_way_compatible


"""
    merge_trees!(ref_tree::T, tree::T)::Tuple{T, Vector{T}} where T<:AbstractNode

Merge two compatible trees, i.e. inserts all cluster of the first tree, which
aren't already in the second tree, into the second tree
"""
function merge_trees!(ref_tree::T, tree::T)::Vector{T} where T<:AbstractNode
    ref_nodes = post_order(ref_tree)
    inserted_nodes = Vector{Node}()
    leaves_dict = Dict{String, Int64}()
    count = 1
    for ref_node in ref_nodes
        if ref_node.nchild == 0
            leaves_dict[ref_node.name] = count
            count += 1
        end # if
    end # for
    nodes = post_order(tree)
    leaves = Vector{Node}()
    cluster_start_indeces = Dict{Node, Int64}()
    for node in nodes
        if node.nchild != 0
            for leaf in leaves
                if (length(node.binary) <= length(leaf.binary) &&
                   node.binary == leaf.binary[1:length(node.binary)])
                    cluster_start_indeces[node] = leaves_dict[leaf.name]
                    break
                end # if
            end # for
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
            push!(leaves, node)
        end # if / else
    end # for
    leaves = order_tree!(tree, cluster_start_indeces)
    set_binary!(tree)
    leaf_ranks = Dict(enumerate([leaf.name for leaf in leaves]))
    leaf_ranks_reverse = Dict(value => key for (key, value) in leaf_ranks)

    xleft_dict, xright_dict = depth_dicts(leaves)

    count = -1
    for ref_node in ref_nodes
        if ref_node.nchild != 0
            start = min_leaf_rank(leaf_ranks_reverse, ref_node)
            stop = max_leaf_rank(leaf_ranks_reverse, ref_node)
            start_node = leaves[start]
            stop_node = leaves[stop]
            xleft = xleft_dict[start_node][1]
            xright = xright_dict[stop_node][1]
            left_path = xleft_dict[start_node][2]
            right_path = xright_dict[stop_node][2]
            if intersect(values(left_path), values(right_path)) == []
                continue
            end
            !xleft.root && delete!(left_path, length(split(xleft.mother.binary, ",")))
            !xright.root && delete!(right_path, length(split(xright.mother.binary,",")))
            depth_left = length(split(xleft.binary, ","))
            depth_right = length(split(xright.binary, ","))
            depth = depth_left >= depth_right
            depth ? p2 = right_path[depth_left] : p2 = left_path[depth_right]
            depth ? p1 = xleft : p1 = xright
            if p2 == p1
                r = p1
                if depth
                    d = left_path[depth_left + 1]
                    e = right_path[depth_left + 1]
                else
                    d = left_path[depth_right + 1]
                    e = right_path[depth_right + 1]
                end # if/else
            elseif p1.mother == p2.mother
                r = p1.mother
                depth ? d = p1 : d = p2
                depth ? e = p2 : e = p1
            else
                continue
            end # if/else
            left = d == r.children[1]
            right = e == r.children[end]
            if !(left && right)
                index_d = findfirst(x -> x == d, r.children)
                index_e = findfirst(x -> x == e, r.children)
                inserted_node =
                    insert_node!(r, r.children[index_d:index_e])
                push!(inserted_nodes, inserted_node)
                # ensures correct depth
                inserted_node.binary = string(r.binary, "z")
                # give unique number to avoid false positive "==" statements
                inserted_node.num = count
                count -= 1
                inserted_depth = length(inserted_node.binary)
                if !left
                    left_path[inserted_depth] = inserted_node
                    right_path[inserted_depth] = inserted_node
                    xleft_dict[start_node] = (inserted_node, left_path)
                end
                if !right
                    right_path[inserted_depth] = inserted_node
                    left_path[inserted_depth] = inserted_node
                    xright_dict[stop_node] = (inserted_node, right_path)
                end
            end # if
        end # if
    end # for
    set_binary!(tree)
    number_nodes!(tree)
    return inserted_nodes
end # function merge_trees!


"""
    order_tree!(root::T, cluster_start_indeces::Dict{T, Int64}, leaves=Vector{T}())
        ::Vector{T} where T<:AbstractNode

Helper function to order a tree based on cluster indeces and return the leaves
of the ordered tree
"""
function order_tree!(root::T, cluster_start_indeces::Dict{T, Int64}, leaves=Vector{T}())::Vector{T} where T<:AbstractNode
    sort!(root.children, by = child -> cluster_start_indeces[child], alg=InsertionSort)
    for child in root.children
        if child.nchild == 0
            push!(leaves, child)
        else
            order_tree!(child, cluster_start_indeces, leaves)
        end # if/else
    end # for
    return leaves
end # function order_tree!


"""
    min_leaf_rank(leaf_ranks::Dict{String, Int64}, node::T)
        ::Int64 where T <: AbstractNode

Recursive helper function to find the lowest ranked leaf descendant of a node
"""
function min_leaf_rank(leaf_ranks::Dict{String, Int64}, node::T)::Int64 where T<:AbstractNode
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
end # function min_leaf_rank


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
end # function max_leaf_rank


"""
    x_left(node::T)::Tuple{T,Vector{T}} where T<:AbstractNode

Helper function to find ancestor of a leaf that has said leaf as leftmost
descendant. Also returns the path from the leaf to the mother of that node.
"""
function x_left(node::T)::Tuple{T, Vector{T}} where T<:AbstractNode
    path = [node]
    while true
        if node.root
            return node, path
        else
            mother = node.mother
            if mother.children[1] != node
                push!(path, node.mother)
                return node, path
            end # if
        node = mother
        push!(path, node)
        end # if/else
    end # while
end # function x_left


"""
    x_right(node::T)::Tuple{T,Vector{T}} where T<:AbstractNode

Helper function to find ancestor of a leaf that has said leaf as rightmost
descendant. Also returns the path from the leaf to the mother of that node.
"""
function x_right(node::T)::Tuple{T, Vector{T}} where T<:AbstractNode
    path = [node]
    while true
        if node.root
            return node, path
        else
            mother = node.mother
            if mother.children[end] != node
                push!(path, node.mother)
                return node, path
            end # if
        node = mother
        push!(path, node)
        end # if/else
    end # while
end # function x_right

function depth_dicts(leaves::Vector{Node})
    xleft_dict = Dict{Node, Tuple{Node, Dict{Int64,Node}}}()
    xright_dict = Dict{Node, Tuple{Node, Dict{Int64, Node}}}()
    for leaf in leaves
        x, path = x_left(leaf)
        xleft_dict[leaf] = (x, Dict(length(split(node.binary, ",")) => node for node in path))
        x, path = x_right(leaf)
        xright_dict[leaf] = (x, Dict(length(split(node.binary, ",")) => node for node in path))
    end
    return xleft_dict, xright_dict
end

"""
    majority_consensus_tree(trees::Vector{T}, percentage::Float64=0.5)::T where T<:AbstractNode

Construct the majority rule consensus tree from a set of trees. By default
includes cluster that occur in over 50% of the trees.
"""
function majority_consensus_tree(trees::Vector{T}, percentage::Float64=0.5)::T where T<:AbstractNode
    first_tree = deepcopy(trees[1])
    leaveset = Set([n.name for n in get_leaves(first_tree)])
    for tree in trees[2:end]
        leaveset2 = Set([n.name for n in get_leaves(tree)])
        leaveset != leaveset2 && throw(ArgumentError("The trees do not have the same set of leaves"))
    end
    nodes = post_order(first_tree)
    leaf_ranks = Dict{String, Int64}()
    count = 0
    # save leaf ranks to order the resulting tree in the end
    for node in nodes
        if node.nchild == 0
            leaf_ranks[node.name] = count
            count += 1
        end # if
    end # for
    ladderize_tree!(first_tree)
    node_counts = convert(Vector{Int64}, ones(length(nodes)))
    count_dict = Dict(zip(nodes, node_counts))
    for tree in trees[2:end]
        nodes = level_order(first_tree)
        is_common_cluster = find_common_clusters(tree, first_tree)
        for node in nodes
            # increment count of clusters of the first tree that are in the other tree
            if is_common_cluster[node.num][1] == true
                count_dict[node] += 1
            # delete clusters which are not
            else
                count_dict[node] -= 1
                if count_dict[node] == 0
                    # delete_node!(first_tree, node)
                    delete_node!(node)
                end # if
            end # else
        end # for
        set_binary!(first_tree)
        number_nodes!(first_tree)
        compatible_tree = one_way_compatible(tree, first_tree)
        inserted_nodes = merge_trees!(compatible_tree, first_tree)
        # intialize counts for the new nodes
        for node in inserted_nodes
            count_dict[node] = 1
        end # for
    end # for
    # find clusters of the final tree in the other tree, store inc_length and count
    nodes = level_order(first_tree)
    count_dict = Dict(zip([n.num for n in nodes], zeros(Int64, length(nodes))))
    inc_length_dict = Dict{Int64, Vector{Float64}}()
    for tree in trees
        is_common_cluster = find_common_clusters(tree, first_tree)
        for node in nodes
            if is_common_cluster[node.num][1]
                count_dict[node.num] += 1
                try
                    push!(inc_length_dict[node.num], is_common_cluster[node.num][2])
                catch KeyError
                    inc_length_dict[node.num] = [is_common_cluster[node.num][2]]
                end # try/catch
            end # if
        end # for
    end # for
    half = length(trees) * percentage
    # delete non-majority clusters
    for node in nodes
        if count_dict[node.num] <= half
            delete_node!(node)
        else
            node.inc_length = mean(inc_length_dict[node.num])
            node.stats["sd"] = std(inc_length_dict[node.num])
            node.stats["median"] = median(inc_length_dict[node.num])
            node.stats["frequency"] = count_dict[node.num] / length(trees)
        end # if
    end # for
    # order the resulting tree
    # """
    nodes = post_order(first_tree)
    cluster_start_indeces = Dict{Node, Int64}()
    for node in nodes
        if node.nchild != 0
            cluster_start_indeces[node] = cluster_start_indeces[node.children[1]]
        else
            cluster_start_indeces[node] = leaf_ranks[node.name]
        end # if/else
    end # for
    order_tree!(first_tree, cluster_start_indeces)
    # """
    set_binary!(first_tree)
    return first_tree
end


"""
    loose_consensus_tree(trees::Vector{T}, percentage::Float64=0.5)::T where T<:AbstractNode

Construct the loose consensus tree from a set of trees. I.e. the clusters appear in at least one
tree and are compatible with all trees.
"""
function loose_consensus_tree(trees::Vector{T}, percentage::Float64=0.5)::T where T<:AbstractNode
    tree_copies = deepcopy(trees)
    r_tree = tree_copies[1]
    leaveset = Set([n.name for n in get_leaves(r_tree)])
    for tree in trees[2:end]
        leaveset2 = Set([n.name for n in get_leaves(tree)])
        leaveset != leaveset2 && throw(ArgumentError("The trees do not have the same set of leaves"))
    end

    nodes = post_order(trees[1])
    leaf_ranks = Dict{String, Int64}()
    count = 0
    # save leaf ranks to order the resulting tree in the end
    for node in nodes
        if node.nchild == 0
            leaf_ranks[node.name] = count
            count += 1
        end # if
    end # for
    for i in 2:length(trees)
        compatible_tree = one_way_compatible(tree_copies[i-1], trees[i])
        merge_trees!(compatible_tree, tree_copies[i])
    end
    loose_tree = tree_copies[end]
    for tree in trees
        loose_tree = one_way_compatible(loose_tree, tree)
    end
    nodes = post_order(loose_tree)
    cluster_start_indeces = Dict{Node, Int64}()
    for node in nodes
        if node.nchild != 0
            cluster_start_indeces[node] = cluster_start_indeces[node.children[1]]
        else
            cluster_start_indeces[node] = leaf_ranks[node.name]
        end # if/else
    end # for
    order_tree!(loose_tree, cluster_start_indeces)
    set_binary!(loose_tree)
    return loose_tree
end

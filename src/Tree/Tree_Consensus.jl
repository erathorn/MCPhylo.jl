"""
    find_common_clusters(ref_tree, tree:T)
        ::Tuple{Vector{Vector{Node}}, Vector{Vector{Node}}} where T<:AbstractNode

Use Day's algorithm to find all clusters of tree2 that also occur in tree1
"""
function find_common_clusters(ref_tree::T, tree::T)::Tuple{Vector{Node}, Vector{Node}} where T<:AbstractNode
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
    occuring_clusters, unique_clusters = Vector{Node}(), Vector{Node}()
    nodes = post_order(tree)
    for node in nodes
        if node.nchild != 0
            leaves = get_leaves(node)
            cluster_indeces = [leaves_dict[leaf.name] for leaf in leaves]
            if length(get_leaves(node)) != maximum(cluster_indeces) -
                                           minimum(cluster_indeces) + 1
                push!(unique_clusters, node)
            elseif haskey(cluster_dict, (minimum(cluster_indeces), maximum(cluster_indeces)))
                push!(occuring_clusters, node)
            else
                push!(unique_clusters, node)
            end # if
        end # if
    end # for
    return occuring_clusters, unique_clusters
end


"""
    are_compatible(cluster_names::Vector{String}, tree::T) where T<:AbstractNode

Check if a cluster of nodes (identified by their names) is compatible with a tree
"""
function are_compatible(cluster::Vector{String}, tree::T)::Bool where T<:AbstractNode
    leaf_names = [leaf.name for leaf in get_leaves(tree)]
    if length(intersect(leaf_names, cluster)) != length(cluster)
        throw(ArgumentError("The cluster contains non-leaf nodes."))
    end # if
    lca = find_lca(tree, cluster)
    children = lca.children
    for child in children
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
    merge_tree(ref_tree::T, tree::T)::T where T<:AbstractNode

Merge two compatible trees
"""
function merge_tree(ref_tree::T, tree::T)::T where T<:AbstractNode
    ref_nodes = post_order(ref_tree)
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
    for node in nodes
        if node.child != 0
            leaves = get_leaves(node)
            cluster_start_indeces[node] = leaves_dict[first(leaves).name]
        else
            cluster_start_indeces[node] = leaves_dict[node.name]
        end # if
    end # for
     leaves = order_tree!(tree, cluster_start_indeces)
     leaf_ranks = Dict(enumerate(leaves))
     leaf_ranks = Dict(values(dict) .=> keys(dict))
end

"""
    order_tree!(root::T, cluster_start_indeces::Dict{T, Int64}, leaves=Vector{String}())
        ::Vector{String} where T<:AbstractNode

Order a tree based on cluster indeces
"""
function order_tree!(root::T, cluster_start_indeces::Dict{T, Int64}, leaves=Vector{String}())::Vector{String} where T<:AbstractNode
    sort!(root.children, by = child -> cluster_start_indeces[child])
    for child in root.children
        if child.nchild == 0
            push!(leaves, child.name)
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
    majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode

Construct the majority rule consensus tree from a set of trees
"""
function majority_consensus_tree(trees::Vector{T})::T where T<:AbstractNode
    all_leaves = [get_leaves(root) for root in trees]
    leaf_names = sort([leaf.name for leaf in all_leaves[1]])
    for leaves in all_leaves[2:end]
        if leaf_names != sort([leaf.name for leaf in leaves])
            throw(ArgumentError("Input Trees need to have identical leaf names"))
        end
    end
    ref_tree = trees[1]
    nodes = level_order(ref_tree)
    node_counts = convert(Vector{Int64}, ones(length(nodes)))
    count_dict = Dict(zip(nodes, node_counts))
    for tree in trees[2:end]
        common_clusters, unique_clusters  = find_common_clusters(tree, ref_tree)
        for node in nodes
            if node in common_clusters
                count_dict[node] += 1
            elseif node.nchild != 0
                count_dict[node] -= 1
                if count_dict[node] == 0
                    remove_node!(ref_tree, node)
                end # if
            end # else
        end # for
    end # for
    return trees[1] # to make it not crash
end

"""
    find_common_clusters(tree:T)
        ::Tuple{Dict{String, Int64}, Dict{Tuple{Int64, Int64}, Int64}} where T<:AbstractNode

Use Day's algorithm to find all clusters of tree2 that also occur in tree1
"""
function find_common_clusters(ref_tree::T, tree::T)::Tuple{Vector{Vector{Node}}, Vector{Vector{Node}}} where T<:AbstractNode
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
    occuring_clusters, missing_clusters = Vector{Vector{Node}}(), Vector{Vector{Node}}()
    nodes = post_order(tree)
    for node in nodes
        if node.nchild != 0
            leaves = get_leaves(node)
            cluster_indeces = [leaves_dict[node.name] for node in leaves]
            if length(get_leaves(node)) != maximum(cluster_indeces) -
                                           minimum(cluster_indeces) + 1
                push!(missing_clusters, leaves)
            elseif haskey(cluster_dict, (minimum(cluster_indeces), maximum(cluster_indeces)))
                push!(occuring_clusters, leaves)
            else
                push!(missing_clusters, leaves)
            end # if
        end # if
    end # for
    return occuring_clusters, missing_clusters
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
        end #if
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
    # TODO: need Day's Algorithm next
    return trees[1] # to make it not crash
end

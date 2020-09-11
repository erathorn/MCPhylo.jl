"""
    SPR(original_root::Node, binary::Bool)::AbstractNode

Performs SPR on tree; takes reference to root of tree, boolean value necessary to determine if tree should be treated as binary or not
Returns reference to root of altered tree
"""
function SPR(original_root::Node,binary::Bool)::AbstractNode
    root = deepcopy(original_root)
    if binary
    spr_tree = perform_spr_binary(root)
else
    spr_tree = perform_spr(root)
end #ifelse
return spr_tree
end
"""
    perform_spr(root::Node)::AbstractNode
performs SPR on non-binary tree

"""
function perform_spr(root::Node)::AbstractNode
    subtree_root, nodes_of_subtree = create_random_subtree(root) #returns reference to subtree to be pruned and reattached
    remove_child!(subtree_root.mother,subtree_root)
    spr_tree = merge_randomly(root,subtree_root) #returns reference to root of tree after reattachment
    return spr_tree
end
"""
    perform_spr_binary(root::Node)::AbstractNode

performs SPR on binary tree

"""
function perform_spr_binary(root::Node)::AbstractNode
    subtree_root, nodes_of_subtree = create_random_subtree(root)#returns reference to subtree to be pruned and reattached
    remove_child!(subtree_root.mother,subtree_root)
    spr_tree = merge_randomly(root,subtree_root)#returns reference to root of tree after reattachment
    return spr_tree
end

"""
    create_random_subtree(root::T) where T<:AbstractNode

selects random, non-root node from tree for use in SPR pruning
"""
function create_random_subtree(root::T)  where T<:AbstractNode
    subtree_root = random_node(root)
    nodes_of_subtree = post_order(subtree_root) #could be used in conjunction with Tree_Pruning.jl
    return subtree_root, nodes_of_subtree
end #function

"""
    merge_randomly(root::T,subtree_root::T)::T where T<:AbstractNode

reattaches subtree to random, non-root node

"""
function merge_randomly(root::T,subtree_root::T)::T  where T<:AbstractNode
    random_mother = random_node(root)
    if random_mother.nchild > 1 #creates "placeholder node" in binary tree if necessary to preserve binarity
        other_child = random_mother.children[2]
        binary_placeholder_node = Node()
        binary_placeholder_node.name = "nameless" #standardizes name of node; constructor defaults to "no_name"

        incoming_length = random_mother.inc_length
        proportion = rand(Uniform(0,1)) #should be a number between 0 and 1
        half_inc_length = incoming_length*proportion
        other_half = incoming_length - half_inc_length #distributes length of node between placeholder node and child to preserve length
        add_child!(binary_placeholder_node,other_child)
        add_child!(binary_placeholder_node,subtree_root)
        random_mother.inc_length = half_inc_length
        binary_placeholder_node.inc_length = other_half

        add_child!(random_mother,binary_placeholder_node)
        remove_child!(random_mother,other_child)
    else
        add_child!(random_mother,subtree_root) #no need for placeholder node if given node has 1 or fewer children
    end #ifelse
    return root
end #function

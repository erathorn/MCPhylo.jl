
function SPR(original_root::Node,binary::Bool)::AbstractNode
    root = deepcopy(original_root)
    if binary
    spr_tree = perform_spr_binary(root)
else
    spr_tree = perform_spr(root)

end #ifelse
return spr_tree
end

function perform_spr(root::Node)::AbstractNode
    subtree_root, nodes_of_subtree = create_random_subtree(root)
    remove_child!(subtree_root.mother,subtree_root)
    spr_tree = merge_randomly(root,subtree_root)
    return spr_tree
end

function perform_spr_binary(root::Node)::AbstractNode
    subtree_root, nodes_of_subtree = create_random_subtree(root)
    remove_child!(subtree_root.mother,subtree_root)
    spr_tree = merge_randomly(root,subtree_root)
    return spr_tree
end


function create_random_subtree(root::T)  where T<:AbstractNode
    subtree_root = random_node(root)
    while subtree_root == root
        subtree_root = random_node(root)
    end #while
    nodes_of_subtree = post_order(subtree_root)
    return subtree_root, nodes_of_subtree
end #function


function merge_randomly(root::T,subtree_root::T)::T  where T<:AbstractNode
    random_mother = random_node(root)
    while random_mother == root
        random_mother = random_node(root)
    end #while
    if random_mother.nchild > 1
        other_child = random_mother.children[2]
        binary_placeholder_node = Node()
        binary_placeholder_node.name = "nameless"

        incoming_length = random_mother.inc_length
        proportion = rand(Uniform(0,1)) #should be a number between 0 and 1
        half_inc_length = incoming_length*proportion
        other_half = incoming_length - half_inc_length
        add_child!(binary_placeholder_node,other_child)
        add_child!(binary_placeholder_node,subtree_root)
        random_mother.inc_length = half_inc_length
        binary_placeholder_node.inc_length = other_half


        add_child!(random_mother,binary_placeholder_node)
        remove_child!(random_mother,other_child)

    else
        add_child!(random_mother,subtree_root)
    end #ifelse
    return root
end #function

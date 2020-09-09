#TODO: BINARY TREES
#TODO: why can not use prune_tree
function SPR(original_root::Node)::AbstractNode
    #TODO: HOW MANY ITERATIONS??
    root = deepcopy(original_root)
    spr_tree = perform_spr(root)
    spr_length = tree_length(spr_tree)
    original_tree_length = tree_length(original_root)

    if spr_length < original_tree_length
        println("GOOD ONE")
        println("Newick Representation of SMALLER TREE")
        println(newick(spr_tree))
        return root
        #perform_spr(spr_tree)
    else
        println("BAD ONE")
        println("before we had this length ", original_tree_length)
        println("and we got this one ",spr_length)
        println("Newick Representation of bigger tree which we got")
        println(newick(spr_tree))
        return root
end #ifelse
    #TODO: IF LENGTH IS SMALLER, THAT'S NEW ORIGINAL TREE
end

function perform_spr(root::Node)::AbstractNode
    subtree_root, nodes_of_subtree = create_random_subtree(root)
    remove_child!(subtree_root.mother,subtree_root)
    spr_tree = merge_randomly(root,subtree_root)
    return spr_tree
end

function test_perform_spr(root::Node)::AbstractNode
    println("BEGINNING TREE OF THIS LOOP: ")
    println(newick(root))
    original_tree_length = tree_length(root)
    subtree_root, nodes_of_subtree = create_random_subtree(root)
    names_of_subtree = [i.name for i in nodes_of_subtree]
    println("THE NAMES TO REMOVE ARE ", names_of_subtree)
    #root_tree_with_no_subtree = prune_tree(root,nodes_of_subtree)
    remove_child!(subtree_root.mother,subtree_root)
    println("PRUNED TREE OF THIS ROOT: ")
    #println(newick(root_tree_with_no_subtree))
    println(newick(root))
    # println("THIS SHOULD NOT HAVE PART")
    # println(newick(root_tree_with_no_subtree))
    spr_tree = merge_randomly(root,subtree_root)
    spr_length = tree_length(spr_tree)
    println("SUBTREE OF THIS LOOP: ")
    println(newick(subtree_root))
    println("RESULT TREE OF THIS LOOP: ")
    println(newick(spr_tree))
    println("ORIGINAL LENGTH: ")
    println(tree_length(root))
    println("RESULT LENGTH: ")
    println(spr_length)
    return spr_tree

end

function create_random_subtree(root::T)  where T<:AbstractNode
    subtree_root = random_node(root)
    #TODO: check how to compare
    while subtree_root == root
        subtree_root = random_node(root)
    end #while
    nodes_of_subtree = post_order(subtree_root)
    return subtree_root, nodes_of_subtree
end #function

function merge_randomly(root::T,subtree_root::T)::T  where T<:AbstractNode
    random_mother = random_node(root)
    add_child!(random_mother,subtree_root)
    return root
end #function

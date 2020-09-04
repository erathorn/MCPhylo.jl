##iterator function:
#input is input tree
#list of nodes used(LoNU) = []
#list of trees(LoT) = [input tree]
#while LoT not empty:
#curtree = LoT.pop
#copyoftree = curtree
#for node x in [all nodes in copyoftree]
#dest = random node from copyoftree (maybe keep track of this somehow 4speed)
#newtree = transformfunction(copyoftree,x, dest)
#copyoftree = curtree
# if newtree.length() < curtree.length()
#LoT.add(newtree)
#
#outside while loop: return curtree
#


##tree transforming function:
#input is root tree, thing to move(TtM), and destination(dest)
# removechild(mother of TtM, TtM)
#IF BINARY:
# make new node, addchild(newnode, TtM)
# addchild(newnode, dest)
#addchild(dest.mother, newnode)
# removechild(dest.mother,dest)
#IF NOT BINARY:
# addchild(dest.mother,TtM)
#

# function transform(root::T, TtM::T, dest::T) where T<:AbstractNode
#     motherofTtM = TtM.mother
#
#TODO: BINARY TREES
#TODO: HOW MANY ITERATIONS??
#TODO: IF LENGTH IS SMALLER, THAT'S NEW ORIGINAL TREE

function SPR(root::Node)::AbstractNode
    println("BEGINNING TREE OF THIS LOOP: ")
    println(newick(root))
    original_tree_length = tree_length(root)
    subtree_root, nodes_of_subtree = create_random_subtree(root)
    names_of_subtree = [i.name for i in nodes_of_subtree]
    println("THE NAMES TO REMOVE ARE ", names_of_subtree)
    root_tree_with_no_subtree = prune_tree(root,names_of_subtree)
    println("PRUNED TREE OF THIS ROOT: ")
    println(newick(root_tree_with_no_subtree))
    # println("THIS SHOULD NOT HAVE PART")
    # println(newick(root_tree_with_no_subtree))
    spr_tree = merge_randomly(root_tree_with_no_subtree,subtree_root)
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
    # if spr_length < original_tree_length
    #     println("WENT TO THE SMALLER TREE")
    #     println("ORIGINAL LENGTH ", original_tree_length)
    #     println("ACHIEVED SMALLER LENGTH ", spr_length)
    #     println()
    #     println()
    #     println()
    #     println("Newick Representation of SMALLER TREE")
    #     println(newick(spr_tree))
    #     SPR(spr_tree)
    # else
    #     ("STARTED AGAIN, LENGTH NOT SMALLER")
    #     println("before we had this length ", original_tree_length)
    #     println("and we got this one ",spr_length)
    #     println("Newick Representation of bigger tree which we got")
    #     println(newick(spr_tree))
    #     SPR(root)

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

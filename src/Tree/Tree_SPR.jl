#include("./Node_Type.jl")
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

function SPR(root::T)::T where T <: AbstractNode
    original_tree_length = tree_length(root)
    subtree_root, nodes_of_subtree = create_random_subtree(root)
    root_tree_with_no_subtree = prune_tree(root,nodes_of_subtree)
    spr_tree = merge_randomly(root_tree_with_no_subtree,subtree_root)
    spr_length = tree_length(spr_tree)
    if spr_length < original_tree_length
        println("CONGRATS")
        println("Original was",original_tree_length)
        println("Now it's",spr_length)
        SPR(spr_tree)
    else
        println("Original was",original_tree_length)
        println("Now it's",spr_length)
        SPR(root)
end # if else
end

function create_random_subtree(root::T)::T  where T<:AbstractNode
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

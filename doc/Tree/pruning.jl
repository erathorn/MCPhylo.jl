include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

tree = MCPhylo.parsing_newick_string("(A,B,(C,D)E)F;")
MCPhylo.number_nodes!(tree)

node_list1 = ["B", "C"]
node_list2 = ["E"]
node_list3 = ["A", "C", "E"]

prune_tree1 = MCPhylo.prune_tree(tree, node_list1)
prune_tree2 = MCPhylo.prune_tree(tree, node_list2)
prune_tree3 = MCPhylo.prune_tree(tree, node_list3)
MCPhylo.prune_tree!(tree, node_list1)

@test newick(prune_tree1) ==  newick(MCPhylo.parsing_newick_string("(A,(D)E)F;"))
@test newick(prune_tree2) ==  newick(MCPhylo.parsing_newick_string("(A,B)F;"))
@test newick(prune_tree3) ==  newick(MCPhylo.parsing_newick_string("(B)F;"))
@test newick(tree) ==  newick(MCPhylo.parsing_newick_string("(A,(D)E)F;"))

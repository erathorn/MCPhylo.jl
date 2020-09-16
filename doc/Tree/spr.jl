using MCPhylo
using Test


binary_tree = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382);")
MCPhylo.number_nodes!(binary_tree)
tree = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);")
MCPhylo.number_nodes!(tree)


spr_binary = MCPhylo.SPR(binary_tree, true)
@test round(MCPhylo.tree_length(binary_tree);digits=3) == round(MCPhylo.tree_length(spr_binary);digits=3)
@test length(MCPhylo.post_order(spr_binary)) == length(MCPhylo.post_order(binary_tree)) ||  length(MCPhylo.post_order(spr_binary)) == length(MCPhylo.post_order(binary_tree)) + 1


spr = MCPhylo.SPR(tree,false)
@test round(MCPhylo.tree_length(tree);digits=3) == round(MCPhylo.tree_length(spr);digits=3)
@test length(MCPhylo.post_order(tree)) == length(MCPhylo.post_order(spr)) || length(MCPhylo.post_order(tree)) + 1 == length(MCPhylo.post_order(spr))

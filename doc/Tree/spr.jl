using MCPhylo
using Test
binary_tree = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382);")
#tree = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);")


spr_binary = MCPhylo.SPR(binary_tree, true)
@test MCPhylo.tree_length(binary_tree) == MCPhylo.tree_length(spr_binary)
@test length(MCPhylo.post_order(binary_tree)) == length(MCPhylo.post_order(binary_tree))

spr = MCPhylo.SPR(tree,false)
@test MCPhylo.tree_length(tree) == MCPhylo.tree_length(spr)
@test length(MCPhylo.post_order(tree)) == length(MCPhylo.post_order(spr))

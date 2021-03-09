using MCPhylo
using Test


binary_tree = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382);")
MCPhylo.number_nodes!(binary_tree)
MCPhylo.set_binary!(binary_tree)

error_tree = MCPhylo.parsing_newick_string("(raccoon:19.19959):0.84600;")
error_tree_not_binary = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300,lizard:5.03):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382);")
tree_binary_two = MCPhylo.parsing_newick_string("((raccoon:19.19959,bear:6.80041)50:0.84600,((sea_lion:11.99700,seal:12.00300)100:7.52973,((monkey:100.85930,cat:47.14069)80:20.59201,weasel:18.87953)75:2.09460)50:3.87382,dog:25.46154);")

@test false == MCPhylo.check_binary(error_tree)
@test false == MCPhylo.check_binary(error_tree_not_binary)
@test_throws ErrorException MCPhylo.SPR(error_tree)
@test true == MCPhylo.check_binary(tree_binary_two)
spr_binary = MCPhylo.SPR(binary_tree)

@test length(Set([n.num for n in post_order(spr_binary)])) == length(Set([n.num for n in post_order(binary_tree)]))

@test round(MCPhylo.tree_length(binary_tree);digits=3) == round(MCPhylo.tree_length(spr_binary);digits=3)

@test length(MCPhylo.post_order(spr_binary)) == length(MCPhylo.post_order(binary_tree))

MCPhylo.number_nodes!(tree_binary_two)

spr_binary_two = MCPhylo.SPR(tree_binary_two)

@test length(Set([n.num for n in post_order(spr_binary_two)])) == length(Set([n.num for n in post_order(tree_binary_two)]))

@test round(MCPhylo.tree_length(tree_binary_two);digits=3) == round(MCPhylo.tree_length(spr_binary_two);digits=3)

@test length(MCPhylo.post_order(spr_binary_two)) == length(MCPhylo.post_order(tree_binary_two))



spr_binary_risky = MCPhylo.risky_SPR(binary_tree)

@test length(Set([n.num for n in post_order(spr_binary_risky)])) == length(Set([n.num for n in post_order(binary_tree)]))
@test round(MCPhylo.tree_length(spr_binary_risky);digits=3) == round(MCPhylo.tree_length(spr_binary);digits=3)
@test length(MCPhylo.post_order(spr_binary_risky)) == length(MCPhylo.post_order(binary_tree))

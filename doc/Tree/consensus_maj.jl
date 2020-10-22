using MCPhylo
using Test

@testset "majority_consensus_tree" begin
    tree1 = MCPhylo.parsing_newick_string("(((A,B),C),(D,E))")
    tree2 = MCPhylo.parsing_newick_string("((A,C),(B,D,E))")
    tree3 = MCPhylo.parsing_newick_string("(((B,C),A),D,E)")
    trees = [tree1, tree2, tree3]
    MCPhylo.number_nodes!.(trees)
    MCPhylo.set_binary!.(trees)
    result = newick(MCPhylo.parsing_newick_string("((A,B,C),D,E)"))
    @test newick(MCPhylo.majority_consensus_tree(trees)) == result
end

@testset "majority_consensus_tree" begin
    tree1 = MCPhylo.parsing_newick_string("(((A,B)AB,C)ABC,(D,E)DE)no_name;")
    tree2 = MCPhylo.parsing_newick_string("((A,C)AC,(B,D,E)BDE)no_name;")
    tree3 = MCPhylo.parsing_newick_string("(((B,C)BC,A)BCA,D,E)no_name;")
    trees = [tree1, tree2, tree3]
    MCPhylo.number_nodes!.(trees)
    MCPhylo.set_binary!.(trees)
    result = newick(MCPhylo.parsing_newick_string("((A,B,C),D,E)"))
    @test newick(MCPhylo.majority_consensus_tree(trees)) == result
end

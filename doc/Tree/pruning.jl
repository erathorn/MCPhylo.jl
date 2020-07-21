include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

tree = MCPhylo.parsing_newick_string("(A,B,(C,D)E)F;")
MCPhylo.number_nodes!(tree)

@testset "Pruning leave"

    node_list2 = ["C"]
    prune_tree2 = MCPhylo.prune_tree(tree, node_list2)

    @test newick(prune_tree2) ==  newick(MCPhylo.parsing_newick_string("(A,B)F;"))


end


@testset "Pruning root"
    node_list1 = ["F"]

    prune_tree1 = MCPhylo.prune_tree(tree, node_list1)

    @test newick(prune_tree1) ==  newick(MCPhylo.parsing_newick_string("(A,(D)E)F;"))
end

@testset "Pruning misc"
    node_list1 = ["B", "C"]
    node_list2 = [ "E"]
    node_list3 = ["A", "C", "E"]

    prune_tree1 = MCPhylo.prune_tree(tree, node_list1)
    prune_tree2 = MCPhylo.prune_tree(tree, node_list2)
    prune_tree3 = MCPhylo.prune_tree(tree, node_list3)


    @test newick(prune_tree1) ==  newick(MCPhylo.parsing_newick_string("(A,(D)E)F;"))
    @test newick(prune_tree2) ==  newick(MCPhylo.parsing_newick_string("(A,B)F;"))
    @test newick(prune_tree3) ==  newick(MCPhylo.parsing_newick_string("(B)F;"))

end


@testset "Pruning inplace"
    node_list1 = ["B", "C"]

    MCPhylo.prune_tree!(tree, node_list1)
    @test newick(tree) ==  newick(MCPhylo.parsing_newick_string("(A,(D)E)F;"))
end

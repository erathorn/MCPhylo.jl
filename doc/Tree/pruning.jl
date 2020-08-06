include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

tree = MCPhylo.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
MCPhylo.number_nodes!(tree)

@testset "Pruning leaf" begin

    node_list = ["C"]
    prune_tree = MCPhylo.prune_tree(tree, node_list)

    @test newick(prune_tree) == newick(MCPhylo.parsing_newick_string("(A,B,((D,E)F)G)H;"))

end

@testset "Pruning root" begin

    node_list = ["H"]

    @test_throws ArgumentError MCPhylo.prune_tree(tree, node_list)

end

@testset "Pruning misc" begin

    node_list1 = ["B", "C"]
    node_list2 = ["F"]
    node_list3 = ["A", "C", "G"]

    prune_tree1 = MCPhylo.prune_tree(tree, node_list1)
    prune_tree2 = MCPhylo.prune_tree(tree, node_list2)
    prune_tree3 = MCPhylo.prune_tree(tree, node_list3)


    @test newick(prune_tree1) == newick(MCPhylo.parsing_newick_string("(A,((D,E)F)G)H;"))
    @test newick(prune_tree2) == newick(MCPhylo.parsing_newick_string("(A,B,(C)G)H;"))
    @test newick(prune_tree3) == newick(MCPhylo.parsing_newick_string("(B)H;"))

end


@testset "Pruning inplace" begin
    node_list1 = ["B", "C"]

    MCPhylo.prune_tree!(tree, node_list1)
    @test newick(tree) ==  newick(MCPhylo.parsing_newick_string("(A,((D,E)F)G)H;"))

end

include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

tree = MCPhylo.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
MCPhylo.number_nodes!(tree)

@testset "find_root" begin
    node1 = find_by_name(tree, "A")
    node2 = find_by_name(tree, "C")
    node3 = find_by_name(tree, "E")
    node4 = find_by_name(tree, "H")

    @test find_root(node1) == tree
    @test find_root(node2) == tree
    @test find_root(node3) == tree
    @test find_root(node4) == tree
end

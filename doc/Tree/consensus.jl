include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

@testset "are_compatible" begin

    tree = MCPhylo.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
    MCPhylo.number_nodes!(tree)
    MCPhylo.set_binary!(tree)

    @testset "are_compatible with nodes" begin
        cluster = [find_by_name(tree, "A")]
        @test MCPhylo.are_compatible(cluster, tree) == true
        cluster = [find_by_name(tree, "A"), find_by_name(tree, "B")]
        @test MCPhylo.are_compatible(cluster, tree) == true
        cluster = [find_by_name(tree, "A"), find_by_name(tree, "C")]
        @test MCPhylo.are_compatible(cluster, tree) == false
        cluster = [find_by_name(tree, "D"), find_by_name(tree, "E")]
        @test MCPhylo.are_compatible(cluster, tree) == true
        cluster = [find_by_name(tree, "C"), find_by_name(tree, "F")]
        @test_throws ArgumentError MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "A"), find_by_name(tree, "B"),
                   find_by_name(tree, "C"), find_by_name(tree, "D"),
                   find_by_name(tree, "E")]
        @test MCPhylo.are_compatible(cluster, tree) == true
        cluster = [find_by_name(tree, "C"), find_by_name(tree, "D"),
                   find_by_name(tree, "E")]
        @test MCPhylo.are_compatible(cluster, tree) == true
    end

    @testset "are_compatible with node names" begin
        cluster = ["C"]
        @test MCPhylo.are_compatible(cluster,tree) == true
        cluster = ["A", "B"]
        @test MCPhylo.are_compatible(cluster,tree) == true
        cluster = ["A", "C"]
        @test MCPhylo.are_compatible(cluster,tree) == false
        cluster = ["D", "E"]
        @test MCPhylo.are_compatible(cluster,tree) == true
        cluster = ["C", "F"]
        @test_throws ArgumentError MCPhylo.are_compatible(cluster,tree)
        cluster = ["A", "B", "C", "D", "E"]
        @test MCPhylo.are_compatible(cluster,tree) == true
        cluster = ["C", "D", "E"]
        @test MCPhylo.are_compatible(cluster,tree) == true
    end
end

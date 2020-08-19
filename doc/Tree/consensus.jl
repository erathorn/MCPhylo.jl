include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

@testset "days_algorithm_preprocess" begin

    ref_tree = MCPhylo.parsing_newick_string("((A,(B,(C,(D,E)))),(F,(G,H)))")
    MCPhylo.number_nodes!(ref_tree)
    MCPhylo.set_binary!(ref_tree)

    tree2 = MCPhylo.parsing_newick_string("((G,(C,(A,(F,E)))),(B,(D,H)))")
    MCPhylo.number_nodes!(tree2)
    MCPhylo.set_binary!(tree2)

    A, B, C, D, E, F, G, H = find_by_name(tree2, "A"), find_by_name(tree2, "B"),
                             find_by_name(tree2, "C"), find_by_name(tree2, "D"),
                             find_by_name(tree2, "E"), find_by_name(tree2, "F"),
                             find_by_name(tree2, "G"), find_by_name(tree2, "H")

    expected_matches = [[G,C,A,F,E,B,D,H]]
    expected_mis_matches = [[F,E], [A,F,E], [C,A,F,E], [G,C,A,F,E], [D,H], [B,D,H]]
    matches, mis_matches = MCPhylo.find_common_clusters(ref_tree, tree2)
    @test matches == expected_matches
    @test mis_matches == expected_mis_matches

    tree3 = MCPhylo.parsing_newick_string("((A,(C,(D,(B,E)))),(G,(F,H)))")
    MCPhylo.number_nodes!(tree3)
    MCPhylo.set_binary!(tree3)

    A, B, C, D, E, F, G, H = find_by_name(tree3, "A"), find_by_name(tree3, "B"),
                             find_by_name(tree3, "C"), find_by_name(tree3, "D"),
                             find_by_name(tree3, "E"), find_by_name(tree3, "F"),
                             find_by_name(tree3, "G"), find_by_name(tree3, "H")

    expected_matches = [[C,D,B,E], [A,C,D,B,E], [G,F,H], [A,C,D,B,E,G,F,H]]
    expected_mis_matches = [[B,E], [D,B,E], [F,H]]
    matches, mis_matches = MCPhylo.find_common_clusters(ref_tree, tree3)
    @test matches == expected_matches
    @test mis_matches == expected_mis_matches

    tree4 = MCPhylo.parsing_newick_string("((A,(B,(C,(D,E)))),(F,(G,H)))")
    MCPhylo.number_nodes!(tree4)
    MCPhylo.set_binary!(tree4)

    A, B, C, D, E, F, G, H = find_by_name(tree4, "A"), find_by_name(tree4, "B"),
                             find_by_name(tree4, "C"), find_by_name(tree4, "D"),
                             find_by_name(tree4, "E"), find_by_name(tree4, "F"),
                             find_by_name(tree4, "G"), find_by_name(tree4, "H")

    expected_matches = [[D,E], [C,D,E], [B,C,D,E], [A,B,C,D,E], [G,H], [F,G,H],
                        [A,B,C,D,E,F,G,H]]
    expected_mis_matches = Vector{Vector{Node}}()
    matches, mis_matches = MCPhylo.find_common_clusters(ref_tree, tree4)
    @test matches == expected_matches
    @test mis_matches == expected_mis_matches

    tree5 = MCPhylo.parsing_newick_string("((G,(X,(A,(F,E)))),(B,(D,H)))")
    MCPhylo.number_nodes!(tree5)
    MCPhylo.set_binary!(tree5)
    @test_throws ArgumentError MCPhylo.find_common_clusters(ref_tree, tree5)

    tree6 = MCPhylo.parsing_newick_string("(X,(G,(C,(A,(F,E)))),(B,(D,H)))))")
    MCPhylo.number_nodes!(tree5)
    MCPhylo.set_binary!(tree5)
    @test_throws ArgumentError MCPhylo.find_common_clusters(ref_tree, tree6)
end


@testset "are_compatible" begin

    tree = MCPhylo.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
    MCPhylo.number_nodes!(tree)
    MCPhylo.set_binary!(tree)

    @testset "are_compatible with nodes" begin
        cluster = [find_by_name(tree, "A")]
        @test MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "A"), find_by_name(tree, "B")]
        @test MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "A"), find_by_name(tree, "C")]
        @test !MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "D"), find_by_name(tree, "E")]
        @test MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "C"), find_by_name(tree, "F")]
        @test_throws ArgumentError MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "A"), find_by_name(tree, "B"),
                   find_by_name(tree, "C"), find_by_name(tree, "D"),
                   find_by_name(tree, "E")]
        @test MCPhylo.are_compatible(cluster, tree)
        cluster = [find_by_name(tree, "C"), find_by_name(tree, "D"),
                   find_by_name(tree, "E")]
        @test MCPhylo.are_compatible(cluster, tree)
    end

    @testset "are_compatible with node names" begin
        cluster = ["C"]
        @test MCPhylo.are_compatible(cluster,tree)
        cluster = ["A", "B"]
        @test MCPhylo.are_compatible(cluster,tree)
        cluster = ["A", "C"]
        @test !MCPhylo.are_compatible(cluster,tree)
        cluster = ["D", "E"]
        @test MCPhylo.are_compatible(cluster,tree)
        cluster = ["C", "F"]
        @test_throws ArgumentError MCPhylo.are_compatible(cluster,tree)
        cluster = ["A", "B", "C", "D", "E"]
        @test MCPhylo.are_compatible(cluster,tree)
        cluster = ["C", "D", "E"]
        @test MCPhylo.are_compatible(cluster,tree)
    end
end


@testset "majority_consensus_tree" begin

    tree1 = MCPhylo.parsing_newick_string("((A,(B,(C,(D,E)))),(F,(G,H)))")
    MCPhylo.number_nodes!(tree1)
    MCPhylo.set_binary!(tree1)
    tree2 = MCPhylo.parsing_newick_string("((G,(C,(A,(F,E)))),(B,(D,H)))")
    MCPhylo.number_nodes!(tree2)
    MCPhylo.set_binary!(tree2)
    tree3 = MCPhylo.parsing_newick_string("((B,(F,(C,(D,G)))),(H,(A,E)))")
    MCPhylo.number_nodes!(tree3)
    MCPhylo.set_binary!(tree3)
    tree4 = MCPhylo.parsing_newick_string("((A,(B,(C,(D,T)))),(F,(G,H)))")
    MCPhylo.number_nodes!(tree4)
    MCPhylo.set_binary!(tree4)

    trees = [tree1, tree2, tree3]
    # TODO: needs to be updated, when majority_consensus_tree method is finished
    @test MCPhylo.majority_consensus_tree(trees) == tree1

    trees = [tree1, tree2, tree3, tree4]
    @test_throws ArgumentError MCPhylo.majority_consensus_tree(trees)

end

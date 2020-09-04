include("../../src/MCPhylo.jl")
using .MCPhylo
using Test

@testset "insert_node" begin
    tree = MCPhylo.parsing_newick_string("(A,B,(C,D)E)F;")
    insert_tree = MCPhylo.parsing_newick_string("(A,B,((C,D)no_name)E)F;")
    insert_tree_newick = newick(insert_tree)
    MCPhylo.number_nodes!(tree)

    children = [find_by_name(tree, "C"), find_by_name(tree, "D")]
    mother = find_by_name(tree, "E")
    insert_node!(mother, children)
    MCPhylo.number_nodes!(tree)
    @test newick(tree) == insert_tree_newick

    MCPhylo.number_nodes!(tree)
    new_node = Node("G")
    append!(children, new_node)
    @test_throws AssertionError insert_node!(mother, children)
 end

@testset "delete_node!" begin
    tree = MCPhylo.parsing_newick_string("(A,B,(C,D)E)F;")
    MCPhylo.number_nodes!(tree)
    remove_tree = MCPhylo.parsing_newick_string("(A,B,C,D)F;")
    remove_tree_newick = newick(remove_tree)
    remove_tree2 = MCPhylo.parsing_newick_string("(B,C,D)F;")
    remove_tree2_newick = newick(remove_tree2)

    delete_node!(find_by_name(tree, "E"))
    @test newick(tree) == remove_tree_newick
    delete_node!(find_by_name(tree, "A"))
    @test newick(tree) == remove_tree2_newick
    @test_throws ArgumentError delete_node!(find_by_name(tree, "F"))
end

@testset "find_lca" begin
    tree = MCPhylo.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
    MCPhylo.set_binary!(tree)

    @testset "find_lca by name" begin
        @test find_lca(tree, ["A","C"]) == tree
        @test find_lca(tree, ["D", "E"]) == find_by_name(tree, "F")
        @test find_lca(tree, ["D", "F"]) == find_by_name(tree, "F")
        @test find_lca(tree, ["C", "D"]) == find_by_name(tree, "G")
        @test find_lca(tree, ["A","B","C"]) == tree
        @test find_lca(tree, ["D", "E", "F", "G"]) == find_by_name(tree, "G")
        @test find_lca(tree, ["D", "E", "C"]) == find_by_name(tree, "G")
    end # find_lca by name

    @testset "find_lca by nodes" begin
        nodes = [find_by_name(tree, "A"), find_by_name(tree, "B"),
                 find_by_name(tree, "C")]
        @test find_lca(tree, nodes) == tree

        nodes = [find_by_name(tree, "D"), find_by_name(tree, "E"),
                 find_by_name(tree, "F"), find_by_name(tree, "G")]
        @test find_lca(tree, nodes) == find_by_name(tree, "G")
    end # find_lca by nodes
end

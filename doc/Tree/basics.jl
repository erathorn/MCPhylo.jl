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

    new_node = Node("G")
    append!(children, new_node)
    @test_throws AssertionError insert_node!(mother, children)

 end

@testset "remove_node" begin

    tree = MCPhylo.parsing_newick_string("(A,B,(C,D)E)F;")
    remove_tree = MCPhylo.parsing_newick_string("(A,B,C,D)F;")
    remove_tree_newick = newick(remove_tree)

    MCPhylo.number_nodes!(tree)

    node_e = find_by_name(tree, "E")
    remove_node!(node_e)
    @test newick(tree) == remove_tree_newick

    node_f = find_by_name(tree, "F")
    @test_throws ArgumentError remove_node!(node_f)

end

@testset "find_lca" begin

    tree = MCPhylo.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")

    @test find_lca(tree, ["A","C"]) == tree
    @test find_lca(tree, ["D", "E"]) == find_by_name(tree, "F")
    @test find_lca(tree, ["D", "F"]) == find_by_name(tree, "G")
    @test find_lca(tree, ["C", "D"]) == find_by_name(tree, "F")

end

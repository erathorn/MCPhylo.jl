using MCPhylo
using Test
filename = "newicktester.nwk"
filepath = joinpath(@__DIR__,filename)
testlist = MCPhylo.load_newick(filepath)
@test testlist == ["(A,B)C;","((A,B)C,(D,E)F)G;"]
testlist = ["(A,B)C;","((A,B)C,(D,E)F)G;"]
for x in testlist
    @test MCPhylo.is_valid_newick_string(x) == true
end
@test MCPhylo.is_valid_newick_string("();") == true

@test MCPhylo.is_valid_newick_string("(A, B)C;") == true
@test MCPhylo.is_valid_newick_string("(A, B )C;") == true
@test MCPhylo.is_valid_newick_string("(A,B));") == false
@test MCPhylo.is_valid_newick_string("invalid") == false

name, len = MCPhylo.parse_name_length("A:0.5")
@test name == "A"
@test len == 0.5
name, len = MCPhylo.parse_name_length("A")
@test name == "A"
@test len == 1.0
name, len = MCPhylo.parse_name_length("")
@test name == "nameless"
@test len == 1.0
name, len = MCPhylo.parse_name_length("A:0.5;")
@test name == "A"
@test len == 0.5

nodes = MCPhylo.parsing_newick_string(testlist[1])
@test nodes.children[1].name == "A"
nodes2 = MCPhylo.parsing_newick_string(testlist[2])
@test nodes2.children[2].name == "F"

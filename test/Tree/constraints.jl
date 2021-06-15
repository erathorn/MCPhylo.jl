using MCPhylo
using Test

@testset "generate_constraints" begin
    dict = generate_constraints(mono=[["A", "B"]], not_mono=[["B", "C"]])
    result1 = Dict([(:mono, [["A","B"]]),
                    (:not_mono, [["B","C"]]),
                    (:exc, Tuple{Vector{String}, Vector{String}}[])])
    @test dict == result1

    generate_constraints!(dict, exc=[(["C", "D"], ["E"])], mono=[["A", "B"]])
    push!(result1[:exc], (["C", "D"], ["E"]))
    @test dict == result1

    dict3 = generate_constraints("./test/Tree/topology.txt")
    dict4 = generate_constraints!(dict, "./test/Tree/topology.txt")
    result3 = Dict([(:mono, [["A","B"]]),
                    (:not_mono, [["C","D"]]),
                    (:exc, [(["E", "F"], ["G"])])])
    result4 = Dict([(:mono, [["A","B"]]),
                    (:not_mono, [["B","C"], ["C", "D"]]),
                    (:exc, [(["C", "D"], ["E"]), (["E", "F"], ["G"])])])
    @test dict3 == result3
    @test dict4 == result4
end


@testset "topological" begin
    tree = MCPhylo.parsing_newick_string("((A,B,(C,D,E)F)G,Z)H;")
end

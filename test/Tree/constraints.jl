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
    result3 = Dict([(:mono, [["A","B"], ["C", "D"]]),
                    (:not_mono, [["C","D", "E"]]),
                    (:exc, [(["E", "F"], ["G"])])])
    result4 = Dict([(:mono, [["A","B"], ["C", "D"]]),
                    (:not_mono, [["B","C"], ["C", "D", "E"]]),
                    (:exc, [(["C", "D"], ["E"]), (["E", "F"], ["G"])])])
    @test dict3 == result3
    @test dict4 == result4

    @test_logs (:warn, "Skipped line with unsupported constraint type 'TEST'.
         Allowed types are 'mono', 'not_mono' and 'exc'")
         generate_constraints("./test/Tree/topology.txt")

    @test_logs (:warn, "Some trivial 'mono' / 'not_mono' type constraints were removed.
         A valid 'mono' / 'not_mono' constraint needs at least 2 elements.")
        generate_constraints(mono=[["A"]])

     @test_logs (:warn, "Some trivial 'exc' type constraints were removed.
      A non-trivial 'exc' constraints needs at least 2 elements in the first, and at least 1 in the second part of the tuple")
        generate_constraints(exc=[(["A"], ["B"])])

    @test_throws FileSyntaxError generate_constraints("./test/Tree/topology2.txt")
    @test_throws FileSyntaxError generate_constraints("./test/Tree/topology3.txt")
end


@testset "topological" begin
    tree = MCPhylo.parsing_newick_string("((A,B,(C,D,E)F)G,(H,I)J)K;")
    constraint1 = generate_constraints(mono=[["A", "B", "C", "D", "E"], ["H", "I"]])
    constraint2 = generate_constraints(mono=[["A", "B"]])
    constraint3 = generate_constraints(not_mono=[["A", "B", "H"], ["C", "D"]])
    constraint4 = generate_constraints(not_mono=[["H", "I"]])
    constraint5 = generate_constraints(exc=[(["A", "B"], ["H"]), (["H", "I"], ["D"])])
    constraint6 = generate_constraints(exc=[(["A", "B"], ["C"])])
    constraint7 = generate_constraints(mono=[["A", "B", "C", "D", "E"], ["H", "I"]],
                                       not_mono=[["A", "B", "H"], ["C", "D"]],
                                       exc=[(["A", "B"], ["H"]), (["H", "I"], ["D"])])
    constraint8 = generate_constraints(mono=[["A", "B"]],
                                       not_mono=[["H", "I"]],
                                       exc=[(["A", "B"], ["C"])])
    @test topological(constraint1, tree)
    @test !(topological(constraint2, tree))
    @test topological(constraint3, tree)
    @test !(topological(constraint4, tree))
    @test topological(constraint5, tree)
    @test !(topological(constraint6, tree))
    @test topological(constraint7, tree)
    @test !(topological(constraint8, tree))
end

@testset "TreeConstraints.jl" begin
    @testset "generate_constraints" begin
        constraints = generate_constraints("./distributions/constraints.txt")
        result = Dict(:exc => [(["E", "F"], ["G"])], 
                      :mono => [["A", "B"], ["C", "D"], ["E", "F"]], 
                      :not_mono => [["C", "D", "E"]])  
        result2 = Dict(:exc => [(["E", "F"], ["G"]), (["A", "E"], ["F"])], 
                       :mono => [["A", "B"], ["C", "D"], ["E", "F"], ["D", "F"]], 
                       :not_mono => [["C", "D", "E"], ["A", "B", "C"]])
        @test constraints[:exc] == result[:exc]
        @test constraints[:mono] == result[:mono]
        @test constraints[:not_mono] == result[:not_mono]
        generate_constraints!(constraints, "./distributions/constraints2.txt")
        @test constraints[:exc] == result2[:exc]
        @test constraints[:mono] == result2[:mono]
        @test constraints[:not_mono] == result2[:not_mono]
        @test_throws FileSyntaxError generate_constraints("./distributions/constraints3.txt")
        @test_throws FileSyntaxError generate_constraints("./distributions/constraints4.txt")
    end # generate_constraints

    @testset "topological" begin

    end
end
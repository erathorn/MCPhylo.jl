@testset "Parser.jl" begin
    @testset "make_tree_with_data" begin

        dim1 = [1.0 1.0 0.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0; 
                0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0]
        dim2 = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0 1.0; 
                0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 1.0 0.0]
        result = make_tree_with_data("./parser/example.csv", "csv", "-", "?", true)
        @test result[2][:,:,1] == dim1
        @test result[2][:,:,2] == dim2
        @test result[1].nchild == 2

        result = make_tree_with_data("./parser/example.csv", "csv", "-", "?", binary=false)
        @test result[1].nchild == 3

        @test_throws FileSyntaxError make_tree_with_data("./parser/example.csv")
        @test_throws ArgumentError make_tree_with_data("./parser/example.csv", "csv")
        @test_throws ArgumentError make_tree_with_data("./parser/example.csv", "csv", "-")
    end

    @testset "ParseNexus == ParseCSV" begin
        @test MCPhylo.ParseNexus("./parser/example2.nex") == MCPhylo.ParseCSV("./parser/example.csv", "-", "?", false)
    end
end
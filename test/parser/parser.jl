@testset "make_tree_with_data" begin
    @test_throws FileSyntaxError make_tree_with_data("test/parser/example.csv")
    @test_throws ArgumentError make_tree_with_data("test/parser/example.csv", "csv")
    @test_throws ArgumentError make_tree_with_data("test/parser/example.csv", "csv", "-")
end
@testset "ParseCSV.jl" begin

    matrix = ['0' '0' '0' '0' '0' '-' '0' '0' '0' '0' '1' '1' '0' '0' '-' '0' '0' '?' '0' '0' '?' '0' '0';
                '0' '-' '0' '0' '0' '0' '?' '0' '0' '0' '0' '?' '?' '0' '0' '?' '0' '0' '0' '1' '?' '?' '0';
                '0' '0' '1' '0' '0' '0' '0' '0' '0' '0' '-' '0' '0' '0' '1' '0' '0' '?' '0' '0' '0' '0' '0']
    langs = ["Swedish_0", "Welsh_N_0", "Sardinian_N_0"]

    @testset "ParseCSV" begin
        data = MCPhylo.ParseCSV("./parser/example.csv", "-", "?") 
        data_target = (3, 23, "-", "?", ["0", "1"], matrix, langs)
        @test data == data_target
    end

    @testset "create_csvdf" begin
        content = ["Swedish_0,0,0,0,0,0,-,0,0,0,0,1,1,0,0,-,0,0,?,0,0,?,0,0", 
                   "Welsh_N_0,0,-,0,0,0,0,?,0,0,0,0,?,?,0,0,?,0,0,0,1,?,?,0", 
                   "Sardinian_N_0,0,0,1,0,0,0,0,0,0,0,-,0,0,0,1,0,0,?,0,0,0,0,0"]
        df = MCPhylo.create_csvdf(content)
        df_target = (langs, matrix)
        @test df == df_target

        # test with different separator
        content2 = ["Swedish_0:0:0:0:0:0:-:0:0:0:0:1:1:0:0:-:0:0:?:0:0:?:0:0",
                   "Welsh_N_0:0:-:0:0:0:0:?:0:0:0:0:?:?:0:0:?:0:0:0:1:?:?:0",
                   "Sardinian_N_0:0:0:1:0:0:0:0:0:0:0:-:0:0:0:1:0:0:?:0:0:0:0:0"]
        df2 = MCPhylo.create_csvdf(content2, ":")
        @test df2 == df_target
    end
end
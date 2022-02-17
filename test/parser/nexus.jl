using Test

@testset "ParseNexus.jl" begin

    matrix = ['1' '1' '1' '0' '0' '0' '0' '0' '0' '-' ;
              '0' '0' '1' '0' '0' '1' '1' '0' '0' '0' ;
              '0' '0' '0' '1' '0' '0' '1' '1' '0' '1' ;
              '1' '0' '0' '-' '-' '0' '1' '0' '0' '1' ;
              '0' '0' '0' '?' '0' '1' '1' '0' '1' '0']
    langs = ["Lang1", "Lang2", "Lang3", "Lang4", "Lang5"]

    @testset "ParseNexus" begin
        nex = MCPhylo.ParseNexus("Example.nex")
        nex_target = (5, 10, "-", "?", ["0", "1"], matrix, langs)
        @test nex == nex_target
        @test_throws FileSyntaxError MCPhylo.ParseNexus("topology.txt")
    end
    
    content = ["DIMENSIONS ntax=5 NCHAR=10;", "FORMAT DATATYPE=Restriction GAP=- MISSING=? interleave=yes;",
                   "MATRIX", "", "Lang1   111000000-", "Lang2   0010011000", "Lang3   0001001101",
                   "Lang4   100--01001", "Lang5   000?011010", ";", "End;"]

    @testset "extract_meta_info" begin
        meta_info = MCPhylo.extract_meta_info(content)
        meta_info_target = (5, 10, "-", "?", "NOSYMBOLS")
        @test meta_info == meta_info_target
    end

    @testset "create_nexusdf" begin
        df = MCPhylo.create_nexusdf(content)
        df_target = (langs, matrix)
        @test df == df_target
    end

    @testset "get_alphabet" begin
        alphabet = MCPhylo.get_alphabet(matrix, "-", "?")
        alphabet_target = ["0", "1"]
        @test alphabet == alphabet_target
    end
end


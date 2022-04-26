@testset "PhyloDist.jl" begin
    
    ntax, nchar, gap, miss, symbols, df, langs = MCPhylo.ParseNexus("./likelihood/primates.nex")

    primates_tree = ParseNewick(
    "(((Tarsius_syrichta:0.0510942,(Lemur_catta:0.0136013,Homo_sapiens:0.0370755)12:0.0343822)
    13:0.224569,(Pan:0.0712342,Gorilla:0.03754)14:0.0295151)15:0.0768634,((Pongo:
    0.020513,Hylobates:0.159117)16:0.239429,Macaca_fuscata:0.454752)
    17:0.0902988,((M_mulatta:0.0644278,M_fascicularis:0.318016)
    18:0.015879,(M_sylvanus:0.100663,Saimiri_sciureus:0.0112774)19:0.2727)20:0.0448203);",
    )

    df = MCPhylo.datafortree(df, langs, primates_tree, symbols, gap, miss)
    pden = ones(4) / 4
    pd = PhyloDist(primates_tree, pden, [1.0], [1.0], JC)

    @testset "PhyloDist" begin
        s = Stochastic(Node(), () -> Normal(0, sqrt(1000)), true)
        s.value = primates_tree

        # Constructors
        pd2 = PhyloDist(s, pden, [1.0], [1.0], JC)
        pd3 = PhyloDist(s, pden, 1.0, 1.0, JC)
        pd4 = PhyloDist(primates_tree, pden, 1.0, 1.0, JC)

        # freeK Constructor
        pd5 = PhyloDist(primates_tree, [1.0], [1.0], freeK)
        @test pd5.substitution_model == freeK
        @test pd5.base_freq == [1.0]
        @test pd5.nbase == 1
        pd6 = PhyloDist(primates_tree, [0.25, 0.25, 0.25, 0.25], [1.0],[1.0], JC)
        # pd5.substitution_model = JC
        # pd5.base_freq = 
        # pd5.nbase = 4

        # all Constructors should produce the same result
        @test all(y -> y == pd, [pd2, pd3, pd4, pd6])

        @test minimum(pd) == -Inf
        @test maximum(pd) == Inf
        @test size(pd) == (4, 1, 22)
    end
    
    @testset "MultiplePhyloDist" begin
        trees = [primates_tree, primates_tree]
        freqs = [0.25 0.25; 0.25 0.25; 0.25 0.25; 0.25 0.25]
        rates = zeros(Float64, 1, 2)
        rates[:] .= 1.0
        
        # Constructors
        mpd = MultiplePhyloDist(trees, freqs, rates, rates, JC)
        mpd2 = MultiplePhyloDist(trees, pden, [1.0], [1.0], JC)
        mpd3 = MultiplePhyloDist(trees, freqs, [1.0], [1.0], JC)
        mpd4 = MultiplePhyloDist(trees, freqs, rates, [1.0], JC)
        mpd5 = MultiplePhyloDist(trees, freqs, [1.0], rates, JC)
        # freeK Constructors
        mpd6 = MultiplePhyloDist(trees, [1.0], [1.0], freeK)
        mpd7 = MultiplePhyloDist(trees, rates, [1.0], freeK)

        # test freeK Constructors
        @test mpd6 == mpd7
        @test mpd6.DistCollector[1].substitution_model == freeK
        @test mpd6.DistCollector[2].substitution_model == freeK
        @test mpd6.DistCollector[1].base_freq == [1.0]
        @test mpd6.DistCollector[2].base_freq == [1.0]
        @test mpd6.DistCollector[1].nbase == 1
        @test mpd6.DistCollector[2].nbase == 1

        # mpd6.DistCollector[1].substitution_model = JC
        # mpd6.DistCollector[2].substitution_model = JC
        # mpd6.DistCollector[1].base_freq = [0.25, 0.25, 0.25, 0.25]
        # mpd6.DistCollector[2].base_freq = [0.25, 0.25, 0.25, 0.25]
        # mpd6.DistCollector[1].nbase = 4
        # mpd6.DistCollector[2].nbase = 4

        # # all Constructors should produce the same result
        # @test all(y -> y == mpd, [mpd2, mpd3, mpd4, mpd5, mpd6])

        # check error throwing
        incompatible_rates = reshape([1.0,2.0,3.0], 1,3)
        @test_throws DimensionMismatch MultiplePhyloDist(trees, [0.25 0.25 0.25 ; 0.25 0.25 0.25], [1.0], [1.0], JC)
        @test_throws DimensionMismatch MultiplePhyloDist(trees, freqs, incompatible_rates, [1.0], JC)
        @test_throws DimensionMismatch MultiplePhyloDist(trees, freqs, [1.0], incompatible_rates, JC)
        @test_throws DimensionMismatch MultiplePhyloDist(trees, incompatible_rates, [1.0], freeK)

        @test minimum(mpd) == -Inf
        @test maximum(mpd) == Inf
        @test size(mpd) == (4, 1, 22, 2)

        # test logpdf functions
        mdf = cat(df, df, dims=4)
        @test logpdf(mpd, mdf) ≈ -8677.360274116634 * 2
        @test MCPhylo.__logpdf(mpd, mdf)[1][1] ≈ -8677.360274116634
        @test MCPhylo.__logpdf(mpd, mdf)[2][1] ≈ -8677.360274116634
    end
end
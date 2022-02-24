@testset "dependent.jl" begin
    s = Stochastic(Node(), () -> Normal(0, sqrt(1000)), true)

    @testset "setmonitor!" begin
        @test_throws BoundsError MCPhylo.setmonitor!(s, [1,2,3])
    end

    @testset "unlist" begin
        @test unlist(s, -11.11) == [-11.11]
    end

    @testset "relist" begin
        @test relist(s, [1,2,3]) == 1
        s = Stochastic(5, () -> Normal(0, sqrt(1000)), true)
        @test relist(s, [1,2,3]) == Array{Int64, 5}(undef, 0, 0, 0, 0, 0)

        l = Logical((x)->x, Node(), true)
        l = Logical(l, 3)
        @test MCPhylo.relistlength(l, [3,2,1]) == (3, 1)
        l = Logical(l, [1,2])
        @test MCPhylo.relistlength(l, [-1, -2, -3]) == ([-1, -2], 2)
        node = Node()
        node.num = 42
        s = Stochastic(Node(), () -> Normal(0, sqrt(1000)), true)
        @test MCPhylo.relistlength(s, node) == (node, 1) 
    end

    @testset "gradlogpdf" begin
        @test gradlogpdf(s, "test") == 0.0
    end

    @testset "Stochastic" begin
        s = Stochastic(() -> Normal(0, sqrt(1000)), Node(), true)
        @test isa(s, Stochastic{GeneralNode{Float64, Int64}})
    end

    @testset "Logical" begin
        l = Logical((x)->x, Node(), true)
        @test isa(l, Logical{GeneralNode{Float64, Int64}})
        @test MCPhylo.ScalarLogical(3.3) == 3.3
        l = Logical(l, 3)
        @test isa(l, Logical{Int64})
        @test l.value == 3
    end

    @testset "logpdf" begin
        # @test logpdf(PNUTSexample[:mtree]) == -25.102234863523872
        # tree = ParseNewick("((A,(B,(C,(D,E)))),(F,(G,H)));")
        # @test MCPhylo.logpdf(PNUTSexample[:mtree], tree) == -35.37127304837715
    end

    @testset "rand" begin
        # result = rand(PNUTSexample[:mypi], 5)
        # @test size(result) == (2, 5) 
        # result = rand(PNUTSexample[:mypi])
        # @test isa(result, Vector{Float64})
        # @test length(result) == 2
    end

    @testset "gradlogpdf" begin
        tree = ParseNewick("((A,(B,(C,(D,E)))),(F,(G,H)));")
        result = (-35.37127304837715, [-0.9571428571428571, -0.9571428571428571, -0.9571428571428571,
                  -0.9571428571428571, -0.9571428571428571, -0.9571428571428571, -0.9571428571428571,
                  -0.9571428571428571, -0.9571428571428571, -0.9571428571428571, -0.9571428571428571, 
                  -0.9571428571428571, -0.9571428571428571, -0.9571428571428571])
        # @test MCPhylo.gradlogpdf(PNUTSexample[:mtree], tree) == result
    end
    @testset "set_inits!" begin
    end



    s = Stochastic(5, () -> Normal(0, sqrt(1000)), true)
    showall(s)
    println()
end
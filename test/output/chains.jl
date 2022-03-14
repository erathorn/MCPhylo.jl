@testset "chains.jl" begin
    c = Chains(reshape(collect(1:60), 10, 2, 3))

    @testset "Chains()" begin
        @test isa(Chains(reshape(collect(1:10), 5, 2)), Chains)
        @test isa(Chains([1.0, 2.0, 3.0]), Chains)
        m = reshape(collect(1:60), 10, 2, 3)
        m2 = Array{String}(undef, 1, 1, 1)
        @test_throws DimensionMismatch Chains(m, m2; names=["C"])
        @test_throws DimensionMismatch Chains(m, m2; chains=[1, 2])
    end

    @testset "getindex" begin
        
        v = c.value
        @test getindex(c, :, :, [true, true, true]).value == v

        @testset "window2inds" begin
            @test dropdims(getindex(c, 1, :, [true, true, true]).value, dims=1) == v[1, :, :]
            @test getindex(c, 4:7, :, [true, true, true]).value == v[4:7, :, :]
            @test getindex(c, 4:2:8, :, [true, true, true]).value == v[4:2:8, :, :]
            @test_throws ArgumentError getindex(c, "Iteration1", :, [true, true, true])
        end

        @testset "names2inds" begin
            @test dropdims(getindex(c, :, 2, [true, true, true]).value, dims=2) == v[:,2,:]
            @test dropdims(getindex(c, :, "Param1", [true, true, true]).value, dims=2) == v[:,1,:]
            @test getindex(c, :, ["Param1", "Param2"], [true, true, true]).value == v
            @test getindex(c, :, [true, true], [true, true, true]).value == v
            @test getindex(c, :, [1, 2], [true, true, true]).value == v
        end

        @testset "iters2inds" begin
            @test MCPhylo.iters2inds(c, "Test") == "Test"
            @test MCPhylo.iters2inds(c, :) == 1:10
            @test MCPhylo.iters2inds(c, 1:10) == 1:1:10
            @test MCPhylo.iters2inds(c, [1.0, 2, 3]) == [1, 2, 3]
        end

        @testset "vcat" begin
            c = Chains(reshape(collect(1:6), 3, 2, 1))
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 4)
            result = vcat(c, c2).value
            target = [1.0 4.0; 2.0 5.0; 3.0 6.0; 7.0 10.0; 8.0 11.0; 9.0 12.0]
            @test dropdims(result, dims=3) == target
            for i in 1:3, j in 1:2
                c.trees[i, j] = "TEST"
                c2.trees[i, j] = "TEST2"
            end
            result = cat(c, c2, dims=1)
            target = ["TEST" "TEST"; "TEST" "TEST"; "TEST" "TEST"; "TEST2" "TEST2"; "TEST2" "TEST2"; "TEST2" "TEST2"]
            @test dropdims(result.trees, dims=3) == target
            c2 = Chains(reshape(collect(7:12), 3, 2, 1))
            @test_throws ArgumentError vcat(c, c2)
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 4, 2)
            @test_throws ArgumentError vcat(c, c2)
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 4, 1, ["C1", "C2"])
            @test_throws ArgumentError vcat(c, c2)
            c2 = Chains(6, 2; start=4, tree_names=[:test])
            @test_throws ArgumentError vcat(c, c2)
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 4, 1, ["Param1", "Param2"], [2])
            @test_throws ArgumentError vcat(c, c2)
        end

        @testset "hcat" begin
            c = Chains(reshape(collect(1:6), 3, 2, 1))
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 1, 1, ["Param3", "Param4"])
            result = hcat(c, c2).value
            target = [1.0 4.0 7.0 10.0; 2.0 5.0 8.0 11.0; 3.0 6.0 9.0 12.0]
            @test dropdims(result, dims=3) == target
            for i in 1:3, j in 1:2
                c.trees[i, j] = "TEST"
                c2.trees[i, j] = "TEST2"
            end
            result = hcat(c, c2)
            target = ["TEST" "TEST" "TEST2" "TEST2"; "TEST" "TEST" "TEST2" "TEST2"; "TEST" "TEST" "TEST2" "TEST2"]
            @test dropdims(result.trees, dims=3) == target
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 4)
            @test_throws ArgumentError hcat(c, c2)
            c2 = Chains(reshape(collect(7:12), 3, 2, 1))
            @test_throws ArgumentError hcat(c, c2)
            c2 = Chains(reshape(collect(7:12), 3, 2, 1), 1, 1, ["Param3", "Param4"], [2])
            @test_throws ArgumentError hcat(c, c2)
        end

        @testset "cat3" begin
            c = Chains(reshape(collect(1:6), 3, 1, 2))
            c2 = Chains(reshape(collect(7:12), 3, 1, 2))
            for i in 1:3, j in 1:2
                c.trees[i, 1, j] = "TEST"
                c2.trees[i, 1, j] = "TEST2"
            end
            result = cat(c, c2, dims=3)
            target = ["TEST" "TEST" "TEST2" "TEST2"; "TEST" "TEST" "TEST2" "TEST2"; "TEST" "TEST" "TEST2" "TEST2"]
            @test dropdims(result.trees, dims=2) == target
        end

        @testset "Base.keys" begin
            c = Chains(reshape(collect(1:6), 3, 2, 1))
            push!(c.names, "TEST")
            @test keys(c) == ["Param1", "Param2", "TEST"]
        end

        @testset "combine" begin
            c = Chains(reshape(collect(1:12), 3, 2, 2))
            result = MCPhylo.combine(c)
            target = [1.0 4.0; 7.0 10.0; 2.0 5.0; 8.0 11.0; 3.0 6.0; 9.0 12.0]
            @test result == target
        end

        @testset "indiscretesupport" begin
            c = Chains(reshape(collect(1:6), 3, 1, 2))
            @test MCPhylo.indiscretesupport(c) == [1]
            c = Chains(reshape(collect(1:0.5:6.5), 3, 2, 2))
            @test MCPhylo.indiscretesupport(c) == [0, 0]
            c = Chains(reshape(collect(1:12), 3, 2, 2))
            @test MCPhylo.indiscretesupport(c, (4, 100)) == [0, 1]
            c = Chains(reshape(collect(1:9), 3, 3, 1))
            @test MCPhylo.indiscretesupport(c, (0, 5)) == [1, 0, 0]
        end

        @testset "link" begin
            c = Chains(reshape(collect(0:0.3:1.5), 2, 3, 1))
            target = [0.0 0.4054651081081642 0.1823215567939546; 0.3 2.1972245773362196 0.4054651081081644]
            @test dropdims(MCPhylo.link(c), dims=3) == target
            c = Chains(reshape(collect(-1:-1:-12), 3, 2, 2))
            @test MCPhylo.link(c) == c.value
        end
    end
end
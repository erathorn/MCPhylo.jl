@testset "chains.jl" begin
    c = Chains(reshape(collect(1:60), 10, 2, 3))
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
    end
end
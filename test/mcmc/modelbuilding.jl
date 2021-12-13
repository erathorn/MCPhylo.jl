
@testset "Logical" begin
    l = Logical(1, (β, xmat) -> xmat * β, false)
    @test isa(l, Logical)
    @test length(l.monitor) == 0
    l = Logical((β, xmat) -> xmat * β, false)
    @test isa(l, Logical)
end

@testset "Stochastic" begin
    s = Stochastic(() -> Normal(0, sqrt(1000)), true)
    @test isa(s, Stochastic)
    @test length(s.monitor) == 1
    s = Stochastic(() -> Normal(0, sqrt(1000)), true)
    @test isa(s, Stochastic)
end
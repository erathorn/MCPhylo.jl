
@testset "Logical" begin
    l = Logical(1, (β, xmat) -> xmat * β, false)
    @test isa(l, Logical)
    @test length(l.monitor) == 0
    l = Logical((β, xmat) -> xmat * β, false)
    @test isa(l, Logical)

    l = Logical((β, xmat) -> xmat * β, 1,  false)
    @test isa(l, Logical)

    s = Logical(Node(), (x)->x, true)
    @test isa(s, Stochastic)
end

@testset "Stochastic" begin
    s = Stochastic(1, () -> Normal(0, sqrt(1000)), true)
    @test isa(s, Stochastic)
    @test length(s.monitor) == 1
    s = Stochastic(() -> Normal(0, sqrt(1000)), true)
    @test isa(s, Stochastic)

    s = Stochastic(() -> Normal(0, sqrt(1000)),1, true)
    @test isa(s, Stochastic)

    s = Stochastic(Node(), ()->UniformBranchLength(), true)
    @test isa(s, Stochastic)
end

@testset"Model" begin
    model = Model(
    y = Stochastic(1, (μ, s2) -> MultivariateNormal(μ, s2), false),
    μ = Logical(1, (β, xmat) -> xmat * β, false),
    β = Stochastic(1, () -> Normal(0, sqrt(1000)), true),
    s2 = Stochastic(() -> InverseGamma(0.001, 0.001), true),
)

    @test isa(mode, Model)
end
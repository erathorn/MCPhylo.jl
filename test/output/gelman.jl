using LinearAlgebra
data = Dict{Symbol,Any}(
    :x => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    :y => [1, 3, 3, 3, 5, 5, 7, 7, 7, 9],
)
data[:xmat] = [ones(10) data[:x]]
model = Model(
    y = Stochastic(1, (μ, s2) -> MultivariateNormal(μ, s2*I), false),
    μ = Logical(1, (β, xmat) -> xmat * β, false),
    β = Stochastic(1, () -> Normal(0, sqrt(1000)), true),
    s2 = Stochastic(() -> InverseGamma(0.001, 0.001), true),
)


@testset "gelmandiag" begin
    Random.seed!(123)
    inits = [
        Dict{Symbol,Any}(
            :y => data[:y],
            :β => rand(Normal(0, 1), 2),
            :s2 => rand(Gamma(1, 1)),
        ) for i = 1:2
    ]

    samplers = [NUTS([:β, :s2], variant = :riemann)]
    setsamplers!(model, samplers)

    sim = mcmc(
        model,
        data,
        inits,
        10000,
        burnin = 500,
        thin = 10,
        chains = 2,
        verbose = false,
    )
    out = gelmandiag(sim).value[1:end-1, :, :]
    minv = minimum(x->isnan(x) ? Inf : x, out)
    maxv = maximum(x->isnan(x) ? -Inf : x, out)
    @test 0.9 < minv < 1.1
    @test 1.0 < maxv < 1.2
end

@testset "heideldiag" begin
    1.0 == heideldiag(randn(500))[2]
end
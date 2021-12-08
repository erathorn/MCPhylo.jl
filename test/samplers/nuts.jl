################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

data = Dict(:x => [1, 2, 3, 4, 5], :y => [1, 3, 3, 3, 5])

model = Model(
    y = Stochastic(1, (μ, s2) -> UnivariateDistribution[Normal(i, s2) for i in μ], false),
    μ = Logical(1, (b0, b1, x) -> b0 .+ b1 .* x, false),
    b0 = Stochastic(() -> Normal(0, 1000), true),
    b1 = Stochastic(() -> Normal(0, 1000), true),
    s2 = Stochastic(() -> InverseGamma(0.001, 0.001), true),
)

inits = [
    Dict{Symbol,Any}(
        :b0 => rand(),
        :b1 => rand(),
        :s2 => rand(),
        :y => data[:y],
        :x => data[:x],
    ) for i = 1:2
]

@testset "Slice" begin
    samplers = [Slice([:b0, :b1, :s2], 1.0)]
    setsamplers!(model, samplers)

    sim = mcmc(model, data, inits, 250, burnin = 125, thin = 1, chains = 2, trees = true)

end

@testset "RWM" begin
    samplers = [RWM([:b0, :b1, :s2], 1.0)]
    setsamplers!(model, samplers)

    sim = mcmc(model, data, inits, 250, burnin = 125, thin = 1, chains = 2, trees = true)

end

@testset "NUTS" begin
    samplers = [NUTS([:b0, :b1, :s2])]
    setsamplers!(model, samplers)

    sim = mcmc(model, data, inits, 250, burnin = 125, thin = 1, chains = 2, trees = true)
end

@testset "HMC" begin
    samplers = [HMC([:b0, :b1, :s2], 0.001, 4)]
    setsamplers!(model, samplers)

    sim = mcmc(model, data, inits, 250, burnin = 125, thin = 1, chains = 2, trees = true)
end

################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

data = Dict{Symbol,Any}(
    :x => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    :y => [1, 3, 3, 3, 5, 5, 7, 7, 7, 9],
)
data[:xmat] = [ones(10) data[:x]]
model = Model(
    y = Stochastic(1, (μ, s2) -> MultivariateNormal(μ, s2), false),
    μ = Logical(1, (β, xmat) -> xmat * β, false),
    β = Stochastic(1, () -> MultivariateNormal(2, sqrt(1000)), true),
    s2 = Stochastic(() -> InverseGamma(0.001, 0.001), true),
)

inits = [
    Dict{Symbol,Any}(:y => data[:y], :β => rand(Normal(0, 1), 2), :s2 => rand(Gamma(1, 1))) for i = 1:2
]

@testset "Slice" begin
    Random.seed!(123)
    samplers = [Slice([:β, :s2], 1.0)]
    setsamplers!(model, samplers)

    sim = mcmc(
        model,
        data,
        inits,
        100000,
        burnin = 50000,
        thin = 10,
        chains = 2,
        verbose=false
    )
    
    r_m = summarystats(sim).value[2:3, 1, 1]
    r_sd = summarystats(sim).value[2:3, 2, 1]
    @test r_m[1]-r_sd[1] <= 0.60 <= r_m[1]+r_sd[1]
    @test r_m[2]-r_sd[2] <= 0.80 <= r_m[2]+r_sd[2]
end

@testset "RWM" begin
    Random.seed!(123)
    samplers = [RWM([:β, :s2], 1.0)]
    setsamplers!(model, samplers)

    sim = mcmc(
        model,
        data,
        inits,
        100000,
        burnin = 50000,
        thin = 10,
        chains = 2,
        verbose=false
    )
    
    r_m = summarystats(sim).value[2:3, 1, 1]
    r_sd = summarystats(sim).value[2:3, 2, 1]
    @test r_m[1]-r_sd[1] <= 0.60 <= r_m[1]+r_sd[1]
    @test r_m[2]-r_sd[2] <= 0.80 <= r_m[2]+r_sd[2]
end

@testset "NUTS_Classic" begin
    Random.seed!(123)
    samplers = [NUTS([:β, :s2], variant=:classic)]
    setsamplers!(model, samplers)

    sim = mcmc(
        model,
        data,
        inits,
        100000,
        burnin = 50000,
        thin = 10,
        chains = 2,
        verbose=false
    )
    
    r_m = summarystats(sim).value[2:3, 1, 1]
    r_sd = summarystats(sim).value[2:3, 2, 1]
    @test r_m[1]-r_sd[1] <= 0.60 <= r_m[1]+r_sd[1]
    @test r_m[2]-r_sd[2] <= 0.80 <= r_m[2]+r_sd[2]
end

@testset "NUTS_Riemann" begin
    Random.seed!(123)
    samplers = [NUTS([:β, :s2], variant=:riemann)]
    setsamplers!(model, samplers)

    sim = mcmc(
        model,
        data,
        inits,
        100000,
        burnin = 50000,
        thin = 10,
        chains = 2,
        verbose=false
    )
    
    r_m = summarystats(sim).value[2:3, 1, 1]
    r_sd = summarystats(sim).value[2:3, 2, 1]
    @test r_m[1]-r_sd[1] <= 0.60 <= r_m[1]+r_sd[1]
    @test r_m[2]-r_sd[2] <= 0.80 <= r_m[2]+r_sd[2]
end

@testset "HMC" begin
    Random.seed!(123)
    samplers = [HMC([:β, :s2], 0.007, 8)]
    setsamplers!(model, samplers)

    sim = mcmc(
        model,
        data,
        inits,
        100000,
        burnin = 50000,
        thin = 10,
        chains = 2,
        verbose=false
    )
    
    r_m = summarystats(sim).value[2:3, 1, 1]
    r_sd = summarystats(sim).value[2:3, 2, 1]
    @test r_m[1]-r_sd[1] <= 0.60 <= r_m[1]+r_sd[1]
    @test r_m[2]-r_sd[2] <= 0.80 <= r_m[2]+r_sd[2]
end

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Serialization


"""
Mamba Julia
Tutorial

"""
data = rand(Normal(0,1), 5000)

my_data=Dict(:data=>data)

model = Model(
     data = Stochastic(1, (μ, σ) -> Normal(μ, σ), false),
        μ = Stochastic(()->Normal(),true),
        σ = Stochastic(()->Exponential(1), true)
)

inits = [Dict(:data => data,
             :μ => randn(),
             :σ => rand()),
        Dict(:data => data,
            :μ => randn(),
            :σ => rand())]

samplers = [NUTS(:μ),
            Slice(:σ, 0.1)]

setsamplers!(model, samplers)

sim = mcmc(model, my_data, inits, 50, burnin=5,thin=5, chains=2)

draw(MCPhylo.plotMC(sim))

gelmandiag(sim)

describe(sim)

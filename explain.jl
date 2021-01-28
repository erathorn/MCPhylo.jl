include("./src/MCPhylo.jl")
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

sim = mcmc(model, my_data, inits, 5000, burnin=500, thin=5, chains=2)

plot(sim, [:contour], var_names=["σ", "likelihood", "μ"], fmt=:pdf, nrow=3, ncol=2)

gelmandiag(sim)

describe(sim)

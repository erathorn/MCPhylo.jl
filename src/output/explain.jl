include("../MCPhylo.jl")
using .MCPhylo
using Serialization

"""
trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

plot1 = Plots.plot(trees[1])
plot2 = Plots.plot(trees[1], treetype=:fan, msc=:blue, mc=:yellow, lc=:white,
           bg=:black, tipfont=(7, :lightgreen))
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

# sim = mcmc(model, my_data, inits, 5000, burnin=500, thin=5, chains=2)
sim = mcmc(model, my_data, inits, 1000, burnin=100, thin=5, chains=2)

plot(sim, [:autocor], fmt=:pdf, nrow=2, ncol=1)



using Revise
using Pkg
Pkg.activate(".")
using MCPhylo
using Random
using LinearAlgebra


Random.seed!(42)

s = rand(InverseWishart(3, diagm(ones(3))))
#s = rand(Exponential())
mu = [0.1,0.4,0.3,0.15,0.05]
μ = randn(3)
#mu1 = mu + rand(5)
#mu1 /= sum(mu1)

my_data = Dict{Symbol, Any}(
  :df => collect(transpose(rand(MultivariateNormal(μ, s), 50)))
);



# model setup
model =  Model(
    df = Stochastic(2, (mu1, cm) ->  MultivariateDistribution[MultivariateNormal(mu1, cm) for i in 1:50], false),
    mu1 = Stochastic(1, () -> Normal(), true),
    cm = Stochastic(2, () -> InverseWishart(3, diagm(ones(3))), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mu1 => randn(3),
    :cm => rand(InverseWishart(3, diagm(ones(3)))),
    :df => my_data[:df],
    
    ) for i in 1:2
    ]

scheme = [#Slice([:mu1, :cm], 1.0),
          #RWM(:mu1, 1.0)
          #SliceSimplex(:mu1, scale = 1)
          NUTS([:mu1, :cm], tree_depth=10, target=0.7)
]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.

sim = mcmc(model, my_data, inits, 250000, burnin=125000,thin=10, chains=2, trees=true)

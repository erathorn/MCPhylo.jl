

using Revise
using Pkg
Pkg.activate(".")
using MCPhylo
using Random
using LinearAlgebra


Random.seed!(42)

s = rand(InverseWishart(50, diagm(ones(50))))
#s = rand(Exponential())
mu = [0.1,0.4,0.3,0.15,0.05]
μ = randn(50)
#mu1 = mu + rand(5)
#mu1 /= sum(mu1)

my_data = Dict{Symbol, Any}(
  :df => collect(transpose(rand(MultivariateNormal(μ, s), 500)))
);



# model setup
model =  Model(
    df = Stochastic(2, (mu1, cm) ->  MultivariateDistribution[MultivariateNormal(mu1, cm) for i in 1:500], false),
    mu1 = Stochastic(1, () -> Normal(), true),
    cm = Stochastic(2, () -> InverseWishart(50, diagm(ones(50))), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mu1 => randn(50),
    :cm => rand(InverseWishart(50, diagm(ones(50)))),
    :df => my_data[:df],
    
    ) for i in 1:2
    ]

scheme = [#Slice([:mu1, :cm], 1.0),
          #RWM(:mu1, 1.0)
          #SliceSimplex(:mu1, scale = 1)
          NUTS([:mu1, :cm], tree_depth=5, target=0.8)
]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.

sim = mcmc(model, my_data, inits, 100, burnin=50,thin=1, chains=1, trees=true)
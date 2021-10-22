

using Revise
using Pkg
Pkg.activate(".")
using MCPhylo
using Random
using LinearAlgebra




s = rand(InverseWishart(10, diagm(ones(5))))
#s = rand(Exponential())
mu = randn(5)

my_data = Dict{Symbol, Any}(
  :df => collect(transpose(rand(MultivariateNormal(mu, s), 500)))
);



# model setup
model =  Model(
    df = Stochastic(2, (mu1, cm) ->  MultivariateDistribution[MultivariateNormal(mu1, cm) for i in 1:500], false, false),
    mu1 = Stochastic(1, () -> Normal(), true),
    cm = Stochastic(2, () -> InverseWishart(10, diagm(ones(5))), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mu1 => randn(5),
    :cm => rand(InverseWishart(10, diagm(ones(5)))),
    :df => my_data[:df],
    
    ) for i in 1:2
    ]

scheme = [#Slice([:mu1, :cm], 1.0)
          RWM([:mu1, :cm], 1.0)
]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.

sim = mcmc(model, my_data, inits, 100, burnin=50,thin=1, chains=1, trees=true)
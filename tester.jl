
#=
tester:
- Julia version: 1.3.1
- Author: erathorn
- Date: 2020-10-07
OPENBLAS_NUM_THREADS=5 JULIA_NUM_THREADS=5 /home/jo/Julia14/julia-1.4.2/bin/julia -O3
=#

using Revise
using Pkg
Pkg.activate(".")
using MCPhylo
using Random
using LinearAlgebra
covmat = rand(InverseWishart(5,Diagonal(ones(5))))
mu = randn(5)
my_data = Dict{Symbol, Any}(
  :df => collect(transpose(rand(MultivariateNormal(mu, covmat), 500))),
  :N => 500
);



# model setup
model =  Model(
    df = Stochastic(2, (mu1, cm, N) ->  MultivariateDistribution[MultivariateNormal(mu1, cm) for i in 1:N], false, false),
    mu1 = Stochastic(1, () -> Normal(), true),
    cm = Stochastic(2, () -> InverseWishart(5,diagm(ones(5))), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mu1 => randn(5),
    :cm => rand(InverseWishart(5,diagm(ones(5)))),
    :df => my_data[:df],
    :N => my_data[:N]
    ) for i in 1:2
    ]

scheme = [MCPhylo.NUTS_Rie([:mu1, :cm], target=0.8, tree_depth=5),
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
#@run 
sim_p = mcmc(model, my_data, inits, 1000, burnin=500,thin=1, chains=2, trees=true)


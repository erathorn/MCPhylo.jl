#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#

extensions = quote

    ## Load needed packages and import methods to be extended

    using Distributions
    import Distributions: minimum, maximum, logpdf

    ## Type declaration
    mutable struct PhyloDist <: ContinuousUnivariateDistribution
        my_tree
        mypi::Real
        data::Array
        #rates::Array
    end
    minimum(d::PhyloDist) = -Inf
    maximum(d::PhyloDist) = Inf

    function logpdf(d::PhyloDist, x::Real)
        rates = ones(3132)

        return MCPhylo.FelsensteinFunction(d.my_tree, d.data, d.mypi, rates,3132)
    end
end
#include("./myMamba.jl")
#using .myMamba
include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo




eval(extensions)
tt, data_arr, df = MCPhylo.make_tree_with_data_mat("./local/IE_Contemporary_Full.nex")


my_data = Dict{Symbol, Any}(
  :mtree => tt,
  :data =>data_arr)




model =  Model(
    y = Stochastic(1,
    (mtree, mypi, data) ->
    begin
        UnivariateDistribution[
        PhyloDist(mtree, mypi, data)]
    end,
    ),
    mypi = Stochastic(
    () -> Truncated(Uniform(0.0,1.0), 0.0, 1.0)
    ),
    mtree = Stochastic(2,
    () -> MCPhylo.CompoundDirichlet(1.0,1.0,0.100,1.0),(false,"false"))

     )
inivals = rand(Uniform(0,1),3132)
inivals =inivals./sum(inivals)

inits = [ Dict(
    :mtree => my_data[:mtree],
    :mypi=> 0.5,
    :y => [-500000],
    :rates=>inivals)]

#scheme = [Slice(:mypi, 0.05), SliceSimplex(:rates)]
scheme = [Slice(:mypi, 0.05),
            #Slice(:blenvec, 0.02),
             MCPhylo.ProbPathHMCSampler(:mtree, 5,5.0)]

setsamplers!(model, scheme)


sim = mcmc(model, my_data, inits, 500, burnin=1, chains=1)



salm = Dict{Symbol, Any}(
  :y => reshape(
    [15, 21, 29, 16, 18, 21, 16, 26, 33, 27, 41, 60, 33, 38, 41, 20, 27, 42],
    3, 6),
  :x => [0, 10, 33, 100, 333, 1000],
  :plate => 3,
  :dose => 6
)


## Model Specification
model = Model(

  y = Stochastic(2,
    (alpha, beta, gamma, x, lambda) ->
      UnivariateDistribution[(
        mu = exp(alpha + beta * log(x[j] + 10) + gamma * x[j] + lambda[i, j]);
        Poisson(mu)) for i in 1:3, j in 1:6
      ],
    false
  ),

  alpha = Stochastic(
    () -> Normal(0, 1000)
  ),

  beta = Stochastic(
    () -> Normal(0, 1000)
  ),

  gamma = Stochastic(
    () -> Normal(0, 1000)
  ),

  lambda = Stochastic(2,
    s2 -> Normal(0, sqrt(s2)),
    false
  ),

  s2 = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  )

)


## Initial Values
inits = [
  Dict(:y => salm[:y], :alpha => 0, :beta => 0, :gamma => 0, :s2 => 10,
       :lambda => zeros(3, 6))
]

## Sampling Scheme
scheme = [Slice([:alpha, :beta, :gamma], [1.0, 1.0, 0.1]),
          AMWG([:lambda, :s2], 0.1)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, salm, inits, 10000, burnin=2500, thin=2, chains=1)
describe(sim)

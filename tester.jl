#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#

extensions = quote

    ## Load needed packages and import methods to be extended
    include("./Likelhood/LikelihoodCalculator.jl")
    using .LikelihoodCalculator
    using Distributions
    import Distributions: minimum, maximum, logpdf, pdf

    ## Type declaration
    mutable struct PhyloDist <: ContinuousUnivariateDistribution
        my_tree::Array
        pi::Real
    end
    minimum(d::PhyloDist) = -Inf
    maximum(d::PhyloDist) = Inf

    function logpdf(d::PhyloDist, x::Real)
        LikelihoodCalculator.FelsensteinFunction(d.my_tree, d.pi)
    end

    
end

include("./Tree/Tree_Module.jl")
include("./Substitution/SubstitutionMat.jl")
include("./Parser/ParseNexus.jl")
using .Tree_Module
using .SubstitutionMat
using .NexusParser
include("./Likelhood/LikelihoodCalculator.jl")
using .LikelihoodCalculator
using DataFrames
using Mamba


eval(extensions)

this_tree = NexusParser.make_tree_with_data("./local/development.nex")

    
my_data = Dict{Symbol, Any}(
  :mtree => Tree_Module.post_order(NexusParser.make_tree_with_data("./local/development.nex")))



model = Model(
    y = Stochastic(1,
    (mtree, pi) -> begin UnivariateDistribution[PhyloDist(mtree, pi)] end
    ),
    pi = Stochastic(
    () -> Truncated(Uniform(0.0,1.0), 0.0, 1.0)
    ) )

inits = [ Dict(:mtree => my_data[:mtree], :pi=> rand(Uniform(0,1)), :y => [-1000000000])]

scheme = [Slice(:pi, 0.05)]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 1000, burnin=200, chains=1)

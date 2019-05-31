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
    mutable struct NewUnivarDist <: ContinuousUnivariateDistribution
        my_tree::Array
        pi::Number
    end
    minimum(d::NewUnivarDist) = -Inf
    maximum(d::NewUnivarDist) = Inf

    function logpdf(d::NewUnivarDist)
        LikelihoodCalculator.FelsensteinFunction(my_tree, pi)
    end

    function pdf(d::NewUnivarDist)
        exp(logpdf(d))
    end
end

include("./Tree/Tree_Module.jl")
include("./Substitution/SubstitutionMat.jl")

include("./Sampler/SamplerFunctions.jl")
include("./Parser/ParseNexus.jl")
using .Tree_Module

using .SubstitutionMat
using .SamplerFunctions
using .NexusParser
using DataFrames
using Mamba
using Distributions

eval(extensions)

this_tree = NexusParser.make_tree_with_data("./local/IE_Contemporary_Full.nex")

data = Dict{Symbol, Any}(
  :mtree => Tree_Module.post_order(this_tree))



model = Model(
    y = Stochastic(1,
    (mtree, pi) -> begin UnivariateDistribution[NewUnivarDist(mtree, pi)] end, false
    ),
    pi = Stochastic(
    () -> Uniform()
    ) )

inits = [ Dict(:mtree => data[:mtree], :pi=> 0.1, :y => [-1000])]

scheme = [Slice(:pi, 0.2)]

setsamplers!(model, scheme)

sim = mcmc(model, data, inits, 10000, burnin=2500)

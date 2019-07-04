#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07
=#

extensions = quote

    ## Load needed packages and import methods to be extended
    include("./Likelihood/LikelihoodCalculator.jl")
    using .LikelihoodCalculator
    using Distributions
    import Distributions: minimum, maximum, logpdf

    ## Type declaration
    mutable struct PhyloDist <: ContinuousUnivariateDistribution
        my_tree::Array
        mypi::Real
        #rates::Array
    end
    minimum(d::PhyloDist) = -Inf
    maximum(d::PhyloDist) = Inf

    function logpdf(d::PhyloDist, x::Real)
        rates = ones(3132)
        return LikelihoodCalculator.FelsensteinFunction(d.my_tree, d.mypi, rates,3132)
    end

end


include("./Tree/Tree_Module.jl")
include("./Substitution/SubstitutionMat.jl")
include("./Parser/ParseNexus.jl")
using .Tree_Module
using .SubstitutionMat
using .NexusParser
include("./Likelihood/LikelihoodCalculator.jl")
using .LikelihoodCalculator
using Mamba


eval(extensions)

this_tree = NexusParser.make_tree_with_data("./local/development.nex")


my_data = Dict{Symbol, Any}(
  :mtree => Tree_Module.post_order(NexusParser.make_tree_with_data("./local/development.nex")))



model = Model(
    y = Stochastic(1,
    (mtree, mypi) ->
    begin UnivariateDistribution[
        PhyloDist(mtree, mypi)]
    end,

    ),
    mypi = Stochastic(
    () -> Truncated(Uniform(0.0,1.0), 0.0, 1.0)
    ),
    #rates = Stochastic(1,
    #()-> Dirichlet(ones(3132))
    #)
     )
inivals = rand(Uniform(0,1),size(this_tree.data)[2])
inivals =inivals./sum(inivals)
inits = [ Dict(:mtree => my_data[:mtree], :mypi=> 0.5, :y => [-500000], :rates=>inivals)]

#scheme = [Slice(:mypi, 0.05), SliceSimplex(:rates)]
scheme = [Slice(:mypi, 0.05)]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 500, burnin=1, chains=1)


x::Array{Float64, (2,2)} = [1. 2.; 1. 2.]

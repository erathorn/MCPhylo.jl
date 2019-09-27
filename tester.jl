#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07

=#

extensions = quote

    ## Load needed packages and import methods to be extended

    using Distributions
    import Distributions: minimum, maximum, logpdf, gradlogpdf

    ## Type declaration
    mutable struct PhyloDist <: ContinuousUnivariateDistribution
        my_tree
        mypi::Real
        data::Array
        rates::Vector
    end
    minimum(d::PhyloDist) = -Inf
    maximum(d::PhyloDist) = Inf

    function logpdf(d::PhyloDist, x::Real)

        mt = MCPhylo.post_order(d.my_tree.value)
        return MCPhylo.FelsensteinFunction(mt, d.mypi, d.rates)
    end

    function gradlogpdf(d::PhyloDist, x::Real)
        mt = MCPhylo.pre_order(d.my_tree.value)
        return MCPhylo.GradiantLog(mt, d.mypi, d.rates)
    end
end

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo

mt = MCPhylo.make_tree_with_data("./local/development.nex")

my_data = Dict{Symbol, Any}(
  :mtree => mt)




model =  Model(
    y = Stochastic(1,
    (mtree, mypi, rates) -> PhyloDist(mtree, mypi, rates)
    ),
    mypi = Stochastic( () -> Truncated(Uniform(0.0,1.0), 0.0, 1.0)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0)),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25])),
     av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]))
     )
inivals = rand(Categorical([0.25, 0.25, 0.25, 0.25]),3132)
inivals2 =rand(Dirichlet([0.25, 0.25, 0.25, 0.25]))

inits = [ Dict(
    :mtree => my_data[:mtree],
    :mypi=> 0.5,
    :y => [-500000],
    :mymap=>inivals,
    :av => inivals2)]


scheme = [ProbPathHMC(:mtree, 3.0,0.2, 0.001, :provided),
         Slice(:mypi, 0.05),
         Slice(:av, 0.02),
         RWM(:mymap, 1)
             ]

setsamplers!(model, scheme)


sim = mcmc(model, my_data, inits, 500, burnin=0, chains=1, trees=true)

MCPhylo.to_file(sim, "")

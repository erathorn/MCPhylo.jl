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
        my_tree::PhyloJul.Node
        mypi::Real
        #rates::Array
    end
    minimum(d::PhyloDist) = -Inf
    maximum(d::PhyloDist) = Inf

    function logpdf(d::PhyloDist, x::Real)
        rates = ones(3132)

        return PhyloJul.FelsensteinFunction(PhyloJul.post_order(d.my_tree), d.mypi, rates,3132)
    end

end
include("PhyloJul.jl")
using .PhyloJul

using Mamba


eval(extensions)

this_tree = PhyloJul.make_tree_with_data("./local/development.nex")


my_data = Dict{Symbol, Any}(
  :mtree => PhyloJul.make_tree_with_data("./local/development.nex"))

my_data[:blenvec] = PhyloJul.get_branchlength_vector!(my_data[:mtree])


model = Model(
    y = Stochastic(1,
    (mtree, blenvec, mypi) ->
    begin

        UnivariateDistribution[
        PhyloDist(PhyloJul.set_branchlength_vector!(mtree, blenvec), mypi)]
    end,

    ),
    mypi = Stochastic(
    () -> Truncated(Uniform(0.0,1.0), 0.0, 1.0)
    ),
    blenvec = Stochastic(1,
        (mtree) -> PhyloJul.CompoundDirichlet(1.0,1.0,0.100,1.0,mtree)


    )#,
    #mtree_po = Logical(
    #(mtree, blenvec) -> Tree_Module.set_branchlength_vector!(mtree, blenvec),
    #false
    #)
    #rates = Stochastic(1,
    #()-> Dirichlet(ones(3132))
    #)
     )
inivals = rand(Uniform(0,1),size(this_tree.data)[2])
inivals =inivals./sum(inivals)

inits = [ Dict(
    :mtree => my_data[:mtree],
    :mypi=> 0.5, :y => [-500000],
    :rates=>inivals,
    :blenvec=>PhyloJul.get_branchlength_vector!(my_data[:mtree]))]

#scheme = [Slice(:mypi, 0.05), SliceSimplex(:rates)]
scheme = [Slice(:mypi, 0.05),
            Slice(:blenvec, 0.02)]

setsamplers!(model, scheme)



sim = mcmc(model, my_data, inits, 500, burnin=1, chains=1)

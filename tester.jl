
#=
tester:
- Julia version: 1.2.0
- Author: erathorn
- Date: 2019-05-07

=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
Random.seed!(1234)
using ForwardDiff

mt, df = make_tree_with_data("local/development.nex") # load your own nexus file



my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :nnodes => size(df)[1],
  :nbase => size(df)[2],
  :nsites => size(df)[3],
)

# model setup
model =  Model(
    df = Stochastic(3,
    (mtree, mypi, rates, nnodes, nbase, nsites) -> PhyloDist(mtree, mypi, rates, nnodes, nbase, nsites), false
    ),
    mypi = Stochastic( () -> Uniform(0.0,1.0)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]))
     )

# intial model values
inivals = rand(Categorical([0.25, 0.25, 0.25, 0.25]),3132)
inivals2 =rand(Dirichlet([1.0, 1.0, 1.0, 1.0]))

inits = [ Dict(
    :mtree => my_data[:mtree],
    :blens => MCPhylo.get_branchlength_vector( my_data[:mtree]),
    :mypi=> 0.5,
    :df => my_data[:df],
    :nnodes => size(my_data[:df])[1],
    :nbase => size(my_data[:df])[2],
    :nsites => size(my_data[:df])[3],
    :mymap=>inivals,
    :av => inivals2
    )
    ]


scheme = [ProbPathHMC(:mtree, 6.0,0.1, 0.001, :provided),
         #BranchSlice(:mtree, 0.05),
         Slice(:mypi, 0.05, Univariate),
         SliceSimplex(:av, scale=0.02),
         RWMC(:mymap)
             ]

setsamplers!(model, scheme)

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 5000, burnin=500,thin=10, chains=1, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "")

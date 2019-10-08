
#=
tester:
- Julia version: 1.2.0
- Author: erathorn
- Date: 2019-05-07

=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo

mt = make_tree_with_data("local/development.nex") # load your own nexus file

my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :nb => length(MCPhylo.get_branchlength_vector(mt)))



# model setup
model =  Model(
    y = Stochastic(1,
    (mtree, mypi, rates) -> PhyloDist(mtree, mypi, rates), false
    ),
    mypi = Stochastic( () -> Uniform(0.0,1.0)),
    #mtreev = Logical(Node(), (mtree, blens) -> MCPhylo.set_branchlength_vector!(mtree, blens)),
    #blens = Stochastic(1, (mtree) -> CompoundDirichletWrap(mtree, 1.0,1.0,0.100,1.0), false),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), false),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]))
     )

# intial model values
inivals = rand(Categorical([0.25, 0.25, 0.25, 0.25]),3132)
inivals2 =rand(Dirichlet([0.25, 0.25, 0.25, 0.25]))

inits = [ Dict(
    :mtree => my_data[:mtree],
    :blens => MCPhylo.get_branchlength_vector( my_data[:mtree]),
    :mypi=> 0.5,
    :y => [-500000],
    :mymap=>inivals,
    :nb => my_data[:nb],
    :av => inivals2)]


scheme = [ProbPathHMC(:mtree, 3.0,0.2, 0.001, :provided),
         BranchSlice(:mtree, 0.05),
         Slice(:mypi, 0.05, Univariate),
         SliceSimplex(:av, scale=0.02),
         RWM(:mymap, 1)
             ]

setsamplers!(model, scheme)

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 5000, burnin=50,thin=10, chains=1, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "")

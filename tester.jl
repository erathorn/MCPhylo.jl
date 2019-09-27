#=
tester:
- Julia version: 1.0.1
- Author: erathorn
- Date: 2019-05-07

=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo

mt = make_tree_with_data("yournexusfile.nex") # load your own nexus file

my_data = Dict{Symbol, Any}(
  :mtree => mt)



# model setup
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

# intial model values
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

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 500, burnin=0, chains=1, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "")

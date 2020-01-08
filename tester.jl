
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



mt, df = make_tree_with_data("local/development.nex"); # load your own nexus file



po = post_order(mt);
for node in po
    node.data = df[:,:,node.num]
end

#mt2 = deepcopy(mt)




my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :nnodes => size(df)[3],
  :nbase => size(df)[1],
  :nsites => size(df)[2],
);

#randomize!(mt2, 1000)
#mt3 = deepcopy(mt)
#mt4 = deepcopy(mt)

# model setup
model =  Model(
    df = Stochastic(3,
    (mtree, mypi, rates, nnodes, nbase, nsites) -> PhyloDist(mtree, mypi, rates, nbase, nsites, nnodes), false, false),
    mypi = Stochastic( () -> Uniform(0,1)),
    mtree = Stochastic(MCPhylo.Node_ncu(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),

    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

# intial model values
#inivals = rand(Categorical([0.25, 0.25, 0.25, 0.25]),3132)
#inivals2 =rand(Dirichlet([1.0, 1.0, 1.0, 1.0]))

inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :mymap=>ones(3132),
    :av => [1,1,1,1]
    ),

    ]


scheme = [PNUTS(:mtree),
          Slice(:mypi, 0.05, Univariate)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
#sim = mcmc(model, my_data, inits, 10000, burnin=1000,thin=50, chains=2, trees=true)
sim = mcmc(model, my_data, inits, 5000, burnin=1000,thin=10, chains=1, trees=true)

sim = mcmc(sim, 10000, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "tneg", 10)

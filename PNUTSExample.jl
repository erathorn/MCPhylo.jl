
#=
tester:
- Julia version: 1.3.1
- Author: erathorn
- Date: 2020-10-07
=#

include("./src/MCPhylo.jl")
using .MCPhylo
using Random

Random.seed!(42)



mt, df = make_tree_with_data("Example.nex", binary=true); # load your own nexus file

po = post_order(mt);
for node in po
    node.data = df[:,:,node.num]
    node.scaler = zeros(1,size(node.data, 2))
end

my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :nnodes => size(df)[3],
  :nbase => size(df)[1],
  :nsites => size(df)[2],
);



# model setup
model =  Model(
    df = Stochastic(3, (mtree, mpi, rates, nnodes, nbase, nsites) ->
                PhyloDist(mtree, mpi, rates, nbase, nsites, nnodes), false, false),
    mpi = Logical( (mypi) -> mypi[1], false),
    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(Dirichlet(2, 1)),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :mymap=>ones(3132),
    :av => [1,1,1,1]
    ),
    ]

scheme = [PNUTS(:mtree),
          SliceSimplex(:mypi)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 500,
    burnin=100,thin=5, chains=1, trees=true)

# request more runs
sim = mcmc(sim, 5000, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "example_run")

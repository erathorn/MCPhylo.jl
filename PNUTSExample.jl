
#=
tester:
- Julia version: 1.3.1
- Author: erathorn
- Date: 2020-10-07
OPENBLAS_NUM_THREADS=5 JULIA_NUM_THREADS=5 /home/jo/Julia14/julia-1.4.2/bin/julia -O3
=#
#using Distributed
#addprocs(2)
using Revise
using Pkg
Pkg.activate(".")
using MCPhylo
using Random
#using LinearAlgebra
#LinearAlgebra.BLAS.set_num_threads(5)

#@everywhere Random.seed!(42)
#
#@everywhere using LinearAlgebra
#@everywhere BLAS.set_num_threads(10)
#include("./src/MCPhylo.jl")
#using .MCPhylo
#using Random
#Random.seed!(42)



#mt, df = make_tree_with_data("Example.nex", binary=true); # load your own nexus file
mt, df = make_tree_with_data("untracked_files/simulation_PNUTS_Paper/out_JC_20-600.nex"); # load your own nexus file



mt2 = deepcopy(mt)
randomize!(mt2)
my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :nnodes => size(df)[3],
  :nbase => size(df)[1],
  :nsites => size(df)[2],
);



# model setup
model =  Model(
    df = Stochastic(3, (mtree, mypi) ->  PhyloDist(mtree, mypi, [1.0], [1.0], JC), false, false),
    mypi = Stochastic(1, () -> Dirichlet(4,1), true),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(Dirichlet(4,1)),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :a => rand(),
    ),
    Dict{Symbol, Union{Any, Real}}(
        :mtree => mt2,
        :mypi=> rand(Dirichlet(4,1)),
        :df => my_data[:df],
        :nnodes => my_data[:nnodes],
        :nbase => my_data[:nbase],
        :nsites => my_data[:nsites],
        :a => rand()
        )
    ]

scheme = [PNUTS(:mtree, target=0.6, targetNNI=4, tree_depth=5, jitter=0.5),
          #PPHMC(:mtree, 0.001, 10, 0.003),
          SliceSimplex(:mypi),
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 10000, burnin=5000,thin=1, chains=2, trees=true)


# request more runs
sim = mcmc(sim, 1000)

# write the output to a path specified as the second argument
to_file(sim, "example_run")


using BenchmarkTools
pd = PhyloDist(mt, [0.25,0.25,0.25,0.25], [1.0], [1.0], JC)

logpdf(pd, df)
@timev logpdf(pd, df)
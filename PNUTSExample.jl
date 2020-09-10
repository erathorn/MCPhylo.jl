
#=
tester:
- Julia version: 1.3.1
- Author: erathorn
- Date: 2020-05-06
=#

include("./src/MCPhylo.jl")
using .MCPhylo
using Random
using Serialization
Random.seed!(42)

mnwk="((Oriya_0:0.002,Italian_0:0.002)14:0.002,((Marathi_0:0.002,Rumanian_List_0:0.002)10:0.002,Welsh_N_0:0.002)13:0.002,(Pashto_0:0.002,(Swedish_0:0.002,(Slovenian_0:0.002,Sardinian_N_0:0.002)11:0.002)12:0.002)15:0.002)16:0.0;"

mt = MCPhylo.parsing_newick_string(mnwk)
MCPhylo.number_nodes!(mt)
MCPhylo.set_binary!(mt)
df = deserialize("untracked_files/devdf.jls")
#mt2, df2 = make_tree_with_data("untracked_files/development.nex", binary=false); # load your own nexus file
#MCPhylo.RF(mt, mt2)
po = post_order(mt);
for node in po
    node.data = log.(df[:,:,node.num])
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
    df = Stochastic(3, (mtree, pi, rates, nnodes, nbase, nsites) -> PhyloDist(mtree, pi, rates, nbase, nsites, nnodes), false, false),
    pi = Logical( (mypi) -> mypi[1], false),
    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> [0.78, 0.22],#rand(Dirichlet(2, 1)),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :mymap=>ones(3132),
    :av => [1,1,1,1]
    ),
    ]

scheme = [#Slice(:mtree, 0.05),
        #    RWM(:mtree, 1),
        PNUTS(:mtree),
          SliceSimplex(:mypi)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 10, burnin=2,thin=2, chains=1, trees=true)

# request more runs
sim = mcmc(sim, 50, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "example_run")

po  = post_order(mt)


blv = get_branchlength_vector(mt)
blv .= 0.002

f(y) = MCPhylo.FelsensteinFunction(po, 0.77, ones(3), df, 3132, y)
f(blv)

using Zygote
r1, r2 = Zygote.pullback(f, blv)
g2 = r2(1.0)[1]
using FiniteDiff
g1 = FiniteDiff.finite_difference_gradient(f, blv)
isapprox(g1, g2)

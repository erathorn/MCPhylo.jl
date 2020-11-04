include("./src/MCPhylo.jl")
using .MCPhylo
using Random


mt_st, df_st = make_tree_with_data("notebook/data-st-64-110.paps.nex"); # load your own nexus file
mt_ie, df_ie = make_tree_with_data("notebook/data-ie-42-208.paps.nex"); # load your own nexus file

function data_to_tree(mt, df)
    po = post_order(mt);
    for node in po
        node.data = df[:,:,node.num]
        node.scaler = zeros(1,size(node.data, 2))
    end
end

data_to_tree(mt_st, df_st)
data_to_tree(mt_ie, df_ie)



my_data = Dict{Symbol, Any}(
  :mtree => [mt_st, mt_ie],
  :df => [df_st, df_ie],
  :nnodes => [size(df_st)[3],size(df_ie)[3]],
  :nbase => [size(df_st)[1],size(df_ie)[1]],
  :nsites => [size(df_st)[2],size(df_ie)[2]]
);



# model setup
model =  Model(
    df_ie = Stochastic(3, (mtree_ie, mypi,  rates, nnodes, nbase, nsites) ->
                            PhyloDist(mtree_ie, mypi[1], ones(my_data[:nsites][2]), my_data[:nbase][2], my_data[:nsites][2], my_data[:nnodes][2]), false, false),
    #df_st = Stochastic(3, (mtree_st, mypi, rates,  nnodes, nbase, nsites) ->
    #                        PhyloDist(mtree_st, mypi[1], ones(my_data[:nsites][1]), my_data[:nbase][1], my_data[:nsites][1], my_data[:nnodes][1]), false, false),
    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree_ie = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes][2]+1, true),
    #mtree_st = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes][1]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree_ie => mt_ie,
    :mtree_st => mt_st,
    :mypi=> rand(Dirichlet(2, 1)),
    :df_ie => my_data[:df][2],
    :df_st => my_data[:df][1],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :mymap=>ones(3132),
    :av => [1,1,1,1]
    ),
    ]

scheme = [RWM(:mtree_ie, :all),
          #PNUTS(:mtree_s t),
          SliceSimplex(:mypi)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 10, burnin=2,thin=1, chains=1, trees=true)

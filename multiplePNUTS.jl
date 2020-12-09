include("./src/MCPhylo.jl")
using .MCPhylo
using Random
using Zygote
using ChainRules
using ChainRulesCore

mt_st, df_st = make_tree_with_data("notebook/data-st-64-110.paps.nex"); # load your own nexus file
mt_ie, df_ie = make_tree_with_data("notebook/data-ie-42-208.paps.nex"); # load your own nexus file
mt_aa, df_aa = make_tree_with_data("notebook/data-aa-58-200.paps.nex"); # load your own nexus file
mt_an, df_an = make_tree_with_data("notebook/data-an-45-210.paps.nex"); # load your own nexus file
mt_pn, df_pn = make_tree_with_data("notebook/data-pn-67-183.paps.nex"); # load your own nexus file

function data_to_tree(mt, df)
    po = post_order(mt);
    for node in po
        node.data = df[:,:,node.num]
        node.scaler = zeros(1,size(node.data, 2))
    end
end

data_to_tree(mt_st, df_st)
data_to_tree(mt_ie, df_ie)
data_to_tree(mt_an, df_an)
data_to_tree(mt_aa, df_aa)
data_to_tree(mt_pn, df_pn)



my_data = Dict{Symbol, Any}(
  :mtree => [mt_st, mt_ie, mt_aa, mt_an, mt_pn],
  :df => [df_st, df_ie, df_aa, df_an, df_pn],
  :nnodes => [size(df_st)[3],size(df_ie)[3],size(df_aa)[3],size(df_an)[3],size(df_pn)[3]],
  :nbase => [size(df_st)[1],size(df_ie)[1],size(df_aa)[1],size(df_an)[1],size(df_pn)[1]],
  :nsites => [size(df_st)[2],size(df_ie)[2],size(df_aa)[2],size(df_an)[2],size(df_pn)[2]],
);


#rates[1:4]
# model setup
model =  Model(
    df_ie = Stochastic(3, (mtree_ie, mypi_ie,  rates) ->
                            PhyloDist(mtree_ie, mypi_ie, rates[1:4], my_data[:nbase][2], my_data[:nsites][2], my_data[:nnodes][2]), false, false),
    #df_st = Stochastic(3, (mtree_st, mypi_st, rates) ->
    #                        PhyloDist(mtree_st, mypi_st, rates[5:8], my_data[:nbase][1], my_data[:nsites][1], my_data[:nnodes][1]), false, false),
    #df_aa = Stochastic(3, (mtree_aa, mypi_aa, rates) ->
    #                        PhyloDist(mtree_aa, mypi_aa, rates[9:12], my_data[:nbase][3], my_data[:nsites][3], my_data[:nnodes][3]), false, false),
    #df_an = Stochastic(3, (mtree_an, mypi_an, rates) ->
    #                        PhyloDist(mtree_an, mypi_an, rates[5:8], my_data[:nbase][4], my_data[:nsites][4], my_data[:nnodes][4]), false, false),
    #df_pn = Stochastic(3, (mtree_pn, mypi_pn, rates) ->
    #                        PhyloDist(mtree_pn, mypi_pn, rates[17:20], my_data[:nbase][5], my_data[:nsites][5], my_data[:nnodes][5]), false, false),
    mypi_ie = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_st = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_aa = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_an = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_pn = Stochastic(1, ()-> Dirichlet(2,1)),
    mtree_ie = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[1],b[1],c[1],d[1]), true),
    #mtree_st = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[2],b[2],c[2],d[2]), true),
    #mtree_aa = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[3],b[3],c[3],d[3]), true),
    #mtree_an = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[4],b[4],c[4],d[4]), true),
    #mtree_pn = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[5],b[5],c[5],d[5]), true),
    a = Stochastic(1, ()->Gamma()),
    b = Stochastic(1, ()->Gamma()),
    c = Stochastic(1, ()->Beta(2,2)),
    d = Stochastic(1, ()->Gamma()),
    rates = Logical(1, (αs) -> vcat(discrete_gamma_rates.(αs, αs, 4)...),false),
    αs = Stochastic(1, () -> Gamma(), true),
     )

# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree_ie => mt_ie,
    :mtree_st => mt_st,
    :mtree_pn => mt_pn,
    :mtree_aa => mt_aa,
    :mtree_an => mt_an,
    :mypi_ie=> rand(Dirichlet(2, 1)),
    :mypi_st=> rand(Dirichlet(2, 1)),
    :mypi_an=> rand(Dirichlet(2, 1)),
    :mypi_aa=> rand(Dirichlet(2, 1)),
    :mypi_pn=> rand(Dirichlet(2, 1)),
    :df_ie => my_data[:df][2],
    :df_st => my_data[:df][1],
    :df_aa => my_data[:df][3],
    :df_an => my_data[:df][4],
    :df_pn => my_data[:df][5],
    :αs => rand(5),
    :co => rand(5),
    :a => rand(5),
    :b => rand(5),
    :c => rand(5),
    :d => rand(5),
    ),
    ]

scheme = [PNUTS(:mtree_ie),
          #PNUTS(:mtree_st),
          #PNUTS(:mtree_an),
          #PNUTS(:mtree_aa),
          #PNUTS(:mtree_pn),
          #RWM(:mtree_ie, :all),
          #RWM(:mtree_st, :all),
          #RWM(:mtree_an, :all),
          #RWM(:mtree_aa, :all),
          #RWM(:mtree_pn, :all),
          #Slice(:mtree_ie, 1.0, Multivariate),
          #Slice(:mtree_an, 1.0, Multivariate),
          SliceSimplex(:mypi_ie),
          #SliceSimplex(:mypi_an),
          #SliceSimplex(:mypi_st),
          #SliceSimplex(:mypi_aa),
          #SliceSimplex(:mypi_pn),
          #Slice(:αs, 1.0, Multivariate),
          #Slice([:a, :b, :c, :d], 1.0, Multivariate)
          #RWM(:αs, 1.0)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 500,
    burnin=100,thin=10, chains=1, trees=true)
to_file(sim, "testMult_")

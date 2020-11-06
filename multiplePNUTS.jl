include("./src/MCPhylo.jl")
using .MCPhylo
using Random


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



# model setup
model =  Model(
    df_ie = Stochastic(3, (mtree_ie, mypi_ie,  rates) ->
                            PhyloDist(mtree_ie, mypi_ie[1], rates[1:4], my_data[:nbase][2], my_data[:nsites][2], my_data[:nnodes][2]), false, false),
    #df_st = Stochastic(3, (mtree_st, mypi_st, rates) ->
    #                        PhyloDist(mtree_st, mypi_st[1], rates[5:8], my_data[:nbase][1], my_data[:nsites][1], my_data[:nnodes][1]), false, false),
    #df_aa = Stochastic(3, (mtree_aa, mypi_aa, rates) ->
    #                        PhyloDist(mtree_aa, mypi_aa[1], rates[9:12], my_data[:nbase][3], my_data[:nsites][3], my_data[:nnodes][3]), false, false),
    #df_an = Stochastic(3, (mtree_an, mypi_an, rates) ->
    #                        PhyloDist(mtree_an, mypi_an[1], rates[13:16], my_data[:nbase][4], my_data[:nsites][4], my_data[:nnodes][4]), false, false),
    #df_pn = Stochastic(3, (mtree_pn, mypi_pn, rates) ->
    #                        PhyloDist(mtree_pn, mypi_pn[1], rates[17:20], my_data[:nbase][5], my_data[:nsites][5], my_data[:nnodes][5]), false, false),
    mypi_ie = Stochastic(1, (co) -> Dirichlet(co)),
    #mypi_st = Stochastic(1, (co) -> Dirichlet(co)),
    #mypi_aa = Stochastic(1, (co) -> Dirichlet(co)),
    #mypi_an = Stochastic(1, (co) -> Dirichlet(co)),
    #mypi_pn = Stochastic(1, (co) -> Dirichlet(co)),
    co = Stochastic(1, () -> Gamma()),
    mtree_ie = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_st = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_aa = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_an = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_pn = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    rates = Logical(1, (αs, βs) -> vcat(discrete_gamma_rates.(αs, βs, 4)...),false),
    αs = Stochastic( () -> Gamma(), true),
    βs = Stochastic( () -> Gamma(), true)
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
    :αs => rand(),
    :βs => rand(),
    :co => rand(2),
    ),
    ]

scheme = [RWM(:mtree_ie, :all),
          #PNUTS(:mtree_st),
          #PNUTS(:mtree_an),
          #PNUTS(:mtree_aa),
          #PNUTS(:mtree_pn),
          #SliceSimplex(:mypi_ie),
          #SliceSimplex(:mypi_an),
          #SliceSimplex(:mypi_st),
          #SliceSimplex(:mypi_aa),
          #SliceSimplex(:mypi_pn),
          #Slice(:co, 1.0, Univariate),
          #Slice(:αs, 1.0, Univariate),
          #Slice(:βs, 1.0, Univariate),
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 10, burnin=5,thin=1, chains=1, trees=true)

to_file(sim, "s")

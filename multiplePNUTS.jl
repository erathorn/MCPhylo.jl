include("./src/MCPhylo.jl")
using .MCPhylo
using Random
Random.seed!(42)

mt_st, df_st = make_tree_with_data("notebook/data-st-64-110.paps.nex"); # load your own nexus file
mt_ie, df_ie = make_tree_with_data("notebook/data-ie-42-208.paps.nex"); # load your own nexus file
mt_aa, df_aa = make_tree_with_data("notebook/data-aa-58-200.paps.nex"); # load your own nexus file
mt_an, df_an = make_tree_with_data("notebook/data-an-45-210.paps.nex"); # load your own nexus file
mt_pn, df_pn = make_tree_with_data("notebook/data-pn-67-183.paps.nex"); # load your own nexus file
mt_st2 = deepcopy(mt_st)
mt_ie2 = deepcopy(mt_ie)
mt_an2 = deepcopy(mt_an)
mt_pn2 = deepcopy(mt_pn)
mt_aa2 = deepcopy(mt_aa)
mt_st3 = deepcopy(mt_st)
mt_ie3 = deepcopy(mt_ie)
mt_aa3 = deepcopy(mt_aa)
mt_an3 = deepcopy(mt_an)
mt_pn3 = deepcopy(mt_pn)

randomize!(mt_st2)
randomize!(mt_ie2)
randomize!(mt_aa2)
randomize!(mt_an2)
randomize!(mt_pn2)

randomize!(mt_st3)
randomize!(mt_ie3)
randomize!(mt_aa3)
randomize!(mt_an3)
randomize!(mt_pn3)


my_data = Dict{Symbol, Any}(
  :mtree => [mt_st, mt_ie, mt_aa, mt_an, mt_pn],
  :df => [df_st, df_ie, df_aa, df_an, df_pn],
  :nnodes => [size(df_st)[3],size(df_ie)[3],size(df_an)[3]],
  :nbase => [size(df_st)[1],size(df_ie)[1],size(df_an)[1]],
  :nsites => [size(df_st)[2],size(df_ie)[2],size(df_an)[2]],
);



# model setup
model =  Model(
    df_ie = Stochastic(3, (tree_ie, pi_ie, rates) -> PhyloDist(tree_ie, pi_ie, [1.0], rates[1:4], Restriction), false, false),
    df_st = Stochastic(3, (tree_st, pi_st, rates) -> PhyloDist(tree_st, pi_st, [1.0], rates[5:8], Restriction), false, false),
    df_aa = Stochastic(3, (tree_aa, pi_aa, rates) -> PhyloDist(tree_aa, pi_aa, [1.0], rates[9:12], Restriction), false, false),
    df_an = Stochastic(3, (tree_an, pi_an, rates) -> PhyloDist(tree_an, pi_an, [1.0], rates[13:16], Restriction), false, false),
    df_pn = Stochastic(3, (tree_pn, pi_pn, rates) -> PhyloDist(tree_pn, pi_pn, [1.0], rates[17:20], Restriction), false, false),
    pi_ie = Stochastic(1, ()-> Dirichlet(2,1)),
    pi_st = Stochastic(1, ()-> Dirichlet(2,1)),
    pi_aa = Stochastic(1, ()-> Dirichlet(2,1)),
    pi_an = Stochastic(1, ()-> Dirichlet(2,1)),
    pi_pn = Stochastic(1, ()-> Dirichlet(2,1)),
    tree_ie = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[1],b[1],c[1],d[1]), true),
    tree_st = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[2],b[2],c[2],d[2]), true),
    tree_aa = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[3],b[3],c[3],d[3]), true),
    tree_an = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[4],b[4],c[4],d[4]), true),
    tree_pn = Stochastic(Node(), (a,b,c,d) -> CompoundDirichlet(a[5],b[5],c[5],d[5]), true),
    a = Stochastic(1, ()->Gamma()),
    b = Stochastic(1, ()->Gamma()),
    c = Stochastic(1, ()->Beta(2,2)),
    d = Stochastic(1, ()->Gamma()),
    rates = Logical(1, (αs) -> vcat(discrete_gamma_rates.(αs, αs, 4)...),false),
    αs = Stochastic(1, () -> Gamma(), true),
     )

# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :tree_ie => mt_ie,
    :tree_st => mt_st,
    :tree_pn => mt_pn,
    :tree_aa => mt_aa,
    :tree_an => mt_an,
    :pi_ie=> rand(Dirichlet(2, 1)),
    :pi_st=> rand(Dirichlet(2, 1)),
    :pi_an=> rand(Dirichlet(2, 1)),
    :pi_aa=> rand(Dirichlet(2, 1)),
    :pi_pn=> rand(Dirichlet(2, 1)),
    :df_ie => my_data[:df][2],
    :df_st => my_data[:df][1],
    :df_aa => my_data[:df][3],
    :df_an => my_data[:df][4],
    :df_pn => my_data[:df][5],
    :αs => rand(5),
    :a => rand(5),
    :b => rand(5),
    :c => rand(5),
    :d => rand(5),
    ),
    Dict{Symbol, Union{Any, Real}}(
    :tree_ie => mt_ie2,
    :tree_st => mt_st2,
    :tree_pn => mt_pn2,
    :tree_aa => mt_aa2,
    :tree_an => mt_an2,
    :pi_ie=> rand(Dirichlet(2, 1)),
    :pi_st=> rand(Dirichlet(2, 1)),
    :pi_an=> rand(Dirichlet(2, 1)),
    :pi_aa=> rand(Dirichlet(2, 1)),
    :pi_pn=> rand(Dirichlet(2, 1)),
    :df_ie => my_data[:df][2],
    :df_st => my_data[:df][1],
    :df_aa => my_data[:df][3],
    :df_an => my_data[:df][4],
    :df_pn => my_data[:df][5],
    :αs => rand(5),
    :a => rand(5),
    :b => rand(5),
    :c => rand(5),
    :d => rand(5),
    ),
    Dict{Symbol, Union{Any, Real}}(
    :tree_ie => mt_ie3,
    :tree_st => mt_st3,
    :tree_pn => mt_pn3,
    :tree_aa => mt_aa3,
    :tree_an => mt_an3,
    :pi_ie=> rand(Dirichlet(2, 1)),
    :pi_st=> rand(Dirichlet(2, 1)),
    :pi_an=> rand(Dirichlet(2, 1)),
    :pi_aa=> rand(Dirichlet(2, 1)),
    :pi_pn=> rand(Dirichlet(2, 1)),
    :df_ie => my_data[:df][2],
    :df_st => my_data[:df][1],
    :df_aa => my_data[:df][3],
    :df_an => my_data[:df][4],
    :df_pn => my_data[:df][5],
    :αs => rand(5),
    :a => rand(5),
    :b => rand(5),
    :c => rand(5),
    :d => rand(5),
    ),
    ]

scheme = [PNUTS(:tree_ie, target=0.8, targetNNI=7),
          PNUTS(:tree_st, target=0.8, targetNNI=7),
          PNUTS(:tree_aa, target=0.8, targetNNI=7),
          PNUTS(:tree_an, target=0.8, targetNNI=7),
          PNUTS(:tree_pn, target=0.8, targetNNI=7),
          SliceSimplex(:pi_ie),
          SliceSimplex(:pi_an),
          SliceSimplex(:pi_st),
          SliceSimplex(:pi_aa),
          SliceSimplex(:pi_pn),
          Slice([:αs, :a, :b, :c, :d], 1.0),
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 250000, burnin=25000,thin=1, chains=3, trees=true)
sim = mcmc(sim, 10, trees=true)
to_file(sim, "testMult_")

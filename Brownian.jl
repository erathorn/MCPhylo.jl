
#=
tester:
- Julia version: 1.3.0
- Author: erathorn
- Date: 2019-05-07
=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using LinearAlgebra
using Random
using StatsFuns
Random.seed!(42)

mt, df = make_tree_with_data("LangData/Dravidian.cc.phy.nex", binary=true); # load your own nexus file
#MCPhylo.rescale_length(mt)

#mvn = MvNormal([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], covariance)
#
#logpdf(mvn, logistic.([1, 1,1,1,0,0,0,0,0]))

l = length(MCPhylo.get_leaves(mt))
n = size(df, 2)
arr = [df[:,:,n.num] for n in MCPhylo.get_leaves(mt)]
arrn = zeros(l, n)
for (ind, i) in enumerate(arr)
    arrn[ind, :] .= i[1,:]
end

po = post_order(mt);
for node in po
    node.data = df[:,:,node.num]
    node.scaler = zeros(1,size(node.data, 2))
end

my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :arrn => arrn,
  :leaves => size(arrn, 1),
  :residuals => size(arrn, 2),
  :nnodes => size(df, 3)
);


model =  Model(
    arrn = Stochastic(2,
    (μi, mtree, P, σii) -> BrownianPhylo(μi, mtree, ones(my_data[:leaves],my_data[:leaves]), σii, P,  my_data[:residuals], my_data[:leaves]),false),
    μi = Logical(2, (μ)->(ones(my_data[:leaves],my_data[:residuals])' .= μ)', false),
    μ = Stochastic(1, (μH, σH) -> Normal(μH, σH), false),
    μH = Stochastic(()->Normal()),
    σH = Stochastic(()->Exponential()),
    P = Stochastic(2, ()  -> Normal(),false),
    #Σ = Logical(2, ζi  -> ζi*ζi',false),
    σii = Logical(1, σi -> σi.^2,false),
    σi = Stochastic(1,(λ) -> Exponential(λ), false),
    #ζi = Logical(1, ζ -> ζ.*my_data[:leaves],false),
    #ζ = Stochastic(1, () -> Dirichlet(my_data[:leaves], 1.0)),
    λ = Stochastic(() -> Exponential()),
    mtree = Stochastic(MCPhylo.Node_ncu(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
     )



# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :arrn => my_data[:arrn],
    :nnodes => my_data[:nnodes],
    :nl => my_data[:nnodes],#my_data[:leaves],
    :nsites => my_data[:residuals],
    :P => randn(my_data[:leaves], my_data[:residuals]),
    :μH => rand(),
    :σH => rand(),
    #:ζ => rand(Dirichlet(my_data[:leaves], 1.0)),
    :Σ => ones(my_data[:leaves],my_data[:leaves]),
    :μ => randn(my_data[:residuals]),
    :λ => rand(),
    :σi => rand(my_data[:residuals]),

    ),
    Dict{Symbol, Union{Any, Real}}(
        :mtree => mt,
        :arrn => my_data[:arrn],
        :nnodes => my_data[:nnodes],
        :nl => my_data[:nnodes],#my_data[:leaves],
        :nsites => my_data[:residuals],
        :P => randn(my_data[:leaves], my_data[:residuals]),
        :μH => rand(),
        :σH => rand(),
        #:ζ => rand(Dirichlet(my_data[:leaves], 1.0)),
        :Σ => ones(my_data[:leaves],my_data[:leaves]),
        :μ => randn(my_data[:residuals]),
        :λ => rand(),
        :σi => rand(my_data[:residuals]),

        )
    ]

scheme = [PNUTS(:mtree),
          #Slice(:mu1, 0.05, Univariate),
          #Slice(:μ, 0.05, Multivariate),

          Slice([:μ, :σi], 0.05, Univariate),
          RWM(:P, 1),
          #SliceSimplex(:ζ),
          #, :λ
          Slice([:σH, :μH, :λ], 0.05, Multivariate)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 10, burnin=2,thin=1, chains=2, trees=true)

#sim = mcmc(sim, 20, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "multiv", 5)

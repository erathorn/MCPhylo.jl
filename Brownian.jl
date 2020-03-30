
#=
tester:
- Julia version: 1.3.0
- Author: erathorn
- Date: 2019-05-07
ohrnstein-uhlenbeck
tree length in compound diriclet
=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using LinearAlgebra
using Random
using StatsFuns
Random.seed!(42)

mt, df = make_tree_with_data("LangData/Dravidian.cc.phy.nex"); # load your own nexus file

#cholesky(MCPhylo.to_covariance(mt).*inits[1][:σi][1]^2)
#inits[1][:σi][1]
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
    (μi, mtree, P, scaler, σii) -> BrownianPhylo(μi, mtree, σii, P, scaler,  my_data[:residuals], my_data[:leaves]),false),

    μi = Logical(2, (μ)->(ones(my_data[:leaves],my_data[:residuals])' .= μ)', false),
    μ = Stochastic(1, (μH, σH) -> Normal(μH, σH), false),

    μH = Stochastic(()->Normal()),
    σH = Stochastic(()->Exponential()),
    P = Stochastic(2, ()  -> Normal(),false),

    σii = Logical(1, σi -> σi.^2,false),
    σi = Stochastic(1,(λ) -> Exponential(λ), false),
    #σi = Stochastic(1,() -> Gamma(2,2), false),
    scaler = Stochastic(()-> Beta(1,1)),
    λ = Stochastic(() -> Exponential()),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),

     )



# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :arrn => my_data[:arrn],
    :nnodes => my_data[:nnodes],
    :nl => my_data[:nnodes],
    :nsites => my_data[:residuals],
    :P => randn(my_data[:leaves], my_data[:residuals]),
    :μH => 1,
    :σH => 1,
    :μ => randn(my_data[:residuals]),
    :λ => 1,
    :σi => rand(my_data[:residuals]),
    :scaler => rand()
    )
    ]

scheme = [#RWM(:mtree, 1),
          #Slice(:mtree, 0.05, Multivariate),
          PNUTS(:mtree),
          #Slice(:μH, 0.05, Univariate),

          #Slice(:scaler, 0.05, Univariate),
          RWM(:P, 1),
          #SliceSimplex(:ζ),
          #, :λ
          Slice(:σi, 0.05, Univariate),
          #:σH,
          Slice([:σH,:μH, :λ, :scaler], 0.05, Multivariate)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 100, burnin=25,thin=5, chains=1, trees=true)
#Juno.@profiler mcmc(model, my_data, inits, 10, burnin=5,thin=1, chains=1, trees=true)
#sim = mcmc(sim, 500, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "multiv", 1)
#sim.trees[end]

#μ = rand(my_data[:residuals])
#μi = (ones(my_data[:leaves],my_data[:residuals])' .= μ)'
#σii = ones(2)
#scaler = 1.0
#P = randn(my_data[:leaves], my_data[:residuals])
#br = BrownianPhylo(μi, mt, σii, P, scaler,  my_data[:residuals], my_data[:leaves])
# #
# blv = get_branchlength_vector(mt)
# blv .*= 0.10
# mt2 = deepcopy(mt)
# set_branchlength_vector!(mt2, blv)
# br2 = BrownianPhylo(μi, mt2, ones(my_data[:leaves],my_data[:leaves]), σii, P,  my_data[:residuals], my_data[:leaves])
#  gradlogpdf(br, arrn)
# r,g=gradlogpdf(br2, arrn)
# # any(isnan.(g))

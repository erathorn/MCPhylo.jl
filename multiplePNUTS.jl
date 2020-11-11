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



# model setup
model =  Model(
    df_ie = Stochastic(3, (mtree_ie, mypi_ie,  rates) ->
                            PhyloDist(mtree_ie, mypi_ie, rates[1:4], my_data[:nbase][2], my_data[:nsites][2], my_data[:nnodes][2]), false, false),
    #df_st = Stochastic(3, (mtree_st, mypi_st, rates) ->
    #                        PhyloDist(mtree_st, mypi_st[1], rates[5:8], my_data[:nbase][1], my_data[:nsites][1], my_data[:nnodes][1]), false, false),
    #df_aa = Stochastic(3, (mtree_aa, mypi_aa, rates) ->
    #                        PhyloDist(mtree_aa, mypi_aa[1], rates[9:12], my_data[:nbase][3], my_data[:nsites][3], my_data[:nnodes][3]), false, false),
    df_an = Stochastic(3, (mtree_an, mypi_an, rates) ->
                            PhyloDist(mtree_an, mypi_an, rates[5:8], my_data[:nbase][4], my_data[:nsites][4], my_data[:nnodes][4]), false, false),
    #df_pn = Stochastic(3, (mtree_pn, mypi_pn, rates) ->
    #                        PhyloDist(mtree_pn, mypi_pn[1], rates[17:20], my_data[:nbase][5], my_data[:nsites][5], my_data[:nnodes][5]), false, false),
    mypi_ie = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_st = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_aa = Stochastic(1, ()-> Dirichlet(2,1)),
    mypi_an = Stochastic(1, ()-> Dirichlet(2,1)),
    #mypi_pn = Stochastic(1, ()-> Dirichlet(2,1)),
    mtree_ie = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_st = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_aa = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    mtree_an = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
    #mtree_pn = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true),
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
    :αs => rand(2),
    :co => rand(2),
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
          Slice(:αs, 1.0, Multivariate)
          #RWM(:αs, 1.0)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 5, burnin=1,thin=1, chains=1, trees=true)


using Cthulhu
using Zygote
using LinearAlgebra
using Juno
using Profile
mr = discrete_gamma_rates(0.5, 0.5, 4)
pd = PhyloDist(mt_ie,[0.5, 0.5], mr, my_data[:nbase][2], my_data[:nsites][2], my_data[:nnodes][2])

#Zygote.@showgrad
@descend mcmc(model, my_data, inits, 5, burnin=1,thin=1, chains=1, trees=true)

r3, r4 = gradlogpdf(pd, df_ie)

r1, r2 = gradlogpdf(pd, df_ie)

blv = get_branchlength_vector(mt_ie)
po = post_order(mt_ie)
Juno.@profiler MCPhylo.FelsensteinFunction(po, [0.5,0.5], mr[1], df_ie, size(df_ie,2), blv)


mu =  1.0 / (2.0 * 0.7 * (1-0.7))

MCPhylo.calc_trans(1.0, 0.7, mu, 0.3)

MA = MCPhylo.calc_trans.(blv, Ref(Array([0.7 0.3])), mu, 0.3)



mr = Zygote._pullback(f, 2.0)
mr2 = Zygote.pullback(g, 2.0)
@code_warntype mr[2](1.0)
mr[2](1.0)
mr2[2](1.0)



function nct(time::Float64, pi_::Array{Float64}, mu, r)
    res = repeat(pi_, 2)
    test = Diagonal(ones(2))
    v = exp(-r*time*mu)
    res = res .*(1-v)
    test = test .* v
    res .+ test
    #res[1,1] = res[1,1]+v
    #res[2,2] = res[2,2]+v
    #res
end

da = rand(2,2157)
da2 = rand(2,2)
lda = [da2, da2, da2]
function test(da, lda)
    res = ones(2,2157)
    for i in 1:3
        res = res .* (lda[i]*da)
    end
    res
end

mtres = mapreduce(y->y*da, (x,z)->x .* z, lda)


using BenchmarkTools
@benchmark Base.maximum(mtres, dims=1)

@benchmark mapslices(maximum, mtres, dims=1)
using LoopVectorization
function mygemmavx!(C, A, B)
   @avx for m ∈ axes(A,1), n ∈ axes(B,2)
       Cmn = zero(eltype(C))
       for k ∈ axes(A,2)
           Cmn += A[m,k] * B[k,n]
       end
       C[m,n] = Cmn
   end
end

mygemmavx!(rand(4,4), rand(4,4), rand(4,4))

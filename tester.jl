
#=
tester:
- Julia version: 1.3.0
- Author: erathorn
- Date: 2019-05-07
=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
Random.seed!(42)

# support for csv
mt, df = make_tree_with_data("LangData/Dravidian.cc.phy.nex"); # load your own nexus file


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
    df = Stochastic(3, (mtree, mypi, rates, nnodes, nbase, nsites) -> PhyloDist(mtree, mypi, rates, nbase, nsites, nnodes), false, false),
    mypi = Stochastic( () -> Uniform(0,1)),
    mtree = Stochastic(MCPhylo.Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )
# number = dimension 0
# vector = dimension 1
# matrix = dimesion 2




# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :mymap=>ones(3132),
    :av => [1,1,1,1]
    ),
    ]

scheme = [PNUTS(:mtree),
          Slice(:mypi, 0.05, Univariate)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
#sim = mcmc(model, my_data, inits, 5000, burnin=1000,thin=5, chains=1, trees=true)

#sim = mcmc(sim, 200, trees=true)

# write the output to a path specified as the second argument
#to_file(sim, "t_dev", 5)
blv = get_branchlength_vector(mt)
using BenchmarkTools
using Zygote

nrec(y) = MCPhylo.FelsensteinFunction(po, pi_, ones(3), df, 231, y)

"""
https://github.com/FluxML/Zygote.jl/issues/292
"""


pi_ = 0.5
nrec(blv)
nrec'(blv)
@benchmark nrec(blv)
@benchmark nrec'($blv)

mt2 = deepcopy(mt)
randomize!(mt2)
po2 = post_order(mt2)
ct = po2

d1 = nrec'(blv)
d2 = nrec'(blv)

function test(mt, pi_, df)
    blv = get_branchlength_vector(mt)
    po = post_order(mt)
    MCPhylo.FelsensteinFunction(po, pi_, ones(3), df, 231, blv)
end

function g_test(mt, pi_, df)
    blv = get_branchlength_vector(mt)
    po = post_order(mt)
    #f(y) = MCPhylo.FelsensteinFunction(po, pi_, ones(3), df, 231, y)

    return Zygote.pullback(MCPhylo.FelsensteinFunction, po, pi_, ones(3), df, 231, blv)
end
Zygote.@code_adjoint MCPhylo.FelsensteinFunction(po, pi_, ones(3), df, 231, blv)
import Zygote.literal_getproperty
Zygote.@adjoint function literal_getproperty( t::S,::Val{T}) where T
	println( "generic literal_getproperty ")
    getproperty(t, T), Δ ->getpropertyAdj(Δ,t,T)
end

ct = mt
ct = mt2
f = g_test(mt, pi_, df)
f2 = g_test(mt2, pi_, df)

mt3 = deepcopy(mt2)
randomize!(mt2)
mytest = f2[2](1.0)[end]
mytest2 = f2[2](1.0)[end]



r = Zygote.pullback(MCPhylo.FelsensteinFunction, po, pi_, ones(3), df, 231, blv)
@benchmark r[2](1.0)[end]
@benchmark g_test(mt, pi_, df)
@benchmark test(mt, pi_, df)


d2 == d1



function rtest(mt, blv, pi_, df)#::Tuple{Float64, Vector{Float64}}
    po = post_order(mt)

    nrec(y) =  MCPhylo.FelsensteinFunction(po, pi_, ones(3), df, 231, y)


    r = Zygote.pullback(nrec, blv)
    #r2 = r[2](1.0)[1]::Vector{Float64}
    #r1 = r[1]::Float64
    r[1],r[2](1.0)
end

rtest(mt, blv, 0.8, df)

pi_ = 0.8


@benchmark rtest($mt, $blv, $pi_, $df)

Zygote.@profile rtest(mt, blv, 0.8, df)

using Juno
@Juno.profiler rtest(mt, blv, 0.8, df)




@benchmark MCPhylo.FelsensteinFunction(mt, 0.8, ones(3), df, 231, blv)
@code_warntype MCPhylo.FelsensteinFunction(po, 0.8, ones(3), df, 231, blv)

r = 1.0
pi_ = 0.8
mu =  1.0 / (2.0 * pi_ * (1-pi_))
mml = MCPhylo.calc_trans.(blv, pi_, mu, r)
od_scaler = zeros(Float64, 1, 231)
@benchmark MCPhylo.Felsenstein_Recursion(mt, 0.8, ones(3), df, 231,mml, od_scaler)


po = post_order(mt)
g(y) = MCPhylo.FelsensteinFunction(po, 0.8, ones(3), df, 231, y)

f(y) = FelsensteinFunction(po, 0.8, ones(3), df, 231, y)
f(blv)
g(blv)

@benchmark g(blv)
@benchmark f(blv)


using Flux
@benchmark Flux.pullback(g, blv)[2](1.0)

using Zygote
@benchmark Zygote.pullback(rec, blv)
@benchmark Zygote.pullback(f, blv)
Zygote.pullback(f, blv)[2](1.0)
z(x,y) = x, y

a = [1,2,3,4]
rest = z.(a, 9)
first.(rest)
last.(rest)


#tree_height(mt)

#for node in MCPhylo.get_leaves(mt)
#    println(node.height)
#end
#unlist(mt)

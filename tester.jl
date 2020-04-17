
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

mt, df = make_tree_with_data("LangData/Dravidian.cc.phy.nex"); # load your own nexus file
#mt, df = make_tree_with_data("local/development.nex"); # load your own nexus file


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
    df = Stochastic(3,
    (mtree, mypi, rates, nnodes, nbase, nsites) -> PhyloDist(mtree, mypi, rates, nbase, nsites, nnodes), false, false),
    mypi = Stochastic( () -> Uniform(0,1)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )



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
sim = mcmc(model, my_data, inits, 5000, burnin=1000,thin=5, chains=1, trees=true)

#sim = mcmc(sim, 1000, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "t_dev", 5)
blv = get_branchlength_vector(mt)
using BenchmarkTools
using Zygote

nrec(y) = MCPhylo.FelsensteinFunction(ct, pi_, ones(3), df, 231, y)

pi_ = 0.5
nrec(blv)
mt2 = deepcopy(mt)
randomize!(mt2)
po2 = post_order(mt2)
ct = po2

d1 = nrec'(blv)
d2 = nrec'(blv)


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


@views function Felsenstein_Recursion(root, pi_::T, rates::Vector{Float64}, data::Array{S,3},
                               n_c::Int64, trmats::Array{Array{S, 2},1},
                               od_scaler::Array{S,2})::Tuple{Array{S,2},Array{S,2}}  where {S<:Real, T<:Real}

    root.nchild == 0 && return data[:,:,root.num], od_scaler

        # here happens the intresting stuff
        child_data::Array{Tuple{Array{S,2},Array{S,2}},1} = Felsenstein_Recursion.(root.children, Ref(pi_), Ref(rates), Ref(data), n_c, Ref(trmats), Ref(od_scaler))

        num_v::Array{Array{S,2},1} = trmats[num_n.(root.children)]
        rv::Array{S,2} = reduce(MCPhylo.pointwise_reduce, num_v .* first.(child_data))
        od_scaler = sum(last.(child_data))
        if !root.root
            scaler = Base.maximum(rv, dims=1)
            od_scaler = od_scaler + log.(scaler)

            rv = rv ./ scaler
        end
        return rv, od_scaler

end


function FelsensteinFunction(root, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Array{S}) where {T<:Real, S<:Real}
    r::Float64 = 1.0
    mu =  1.0 / (2.0 * pi_ * (1-pi_))
    mml = MCPhylo.calc_trans.(blv, pi_, mu, r)
    od_scaler = zeros(S, 1, n_c)

    rv, od_scaler = Felsenstein_Recursion(root, pi_, rates, data, n_c, mml, od_scaler)


    sum(log.(sum(rv .* Array([pi_, 1.0-pi_]), dims=1)) + od_scaler)
end

simi

dataf(node)=node.data


num_n(node) = node.num
function FelsensteinFunction(tree_postorder::Vector, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Array{S}) where {S<:Real, T<:Real}
    r::Float64 = 1.0
    mu =  1.0 / (2.0 * pi_ * (1-pi_))
    mml = MCPhylo.calc_trans.(blv, pi_, mu, r)

    root_node = last(tree_postorder)
    #root_node.scaler = zeros(size(root_node.scaler))

    #datac = zero(data)
    rns = zeros(S, 1, n_c)
    @views @inbounds for node in tree_postorder
        if node.nchild > 0
            res = node_loop(node, mml)
            if !node.root
                scaler = Base.maximum(res, dims=1)
                rns += log.(scaler)
                node.data = res ./ scaler
            end #if
        end #if
    end # for

    return likelihood_root(root_node, pi_, rns)#res
end # function

@inline function likelihood_root(root, pi_::S, rns::Array{S,2})::Real where {S<:Number}
    sum(log.(sum(root.data::Array{S,2} .* Array([pi_, 1.0-pi_]), dims=1))+rns)
end

function node_loop(node::T, mml::Array{Array{S, 2},1})::Array{S,2} where {T, S<:Real}
    reduce(pointwise_reduce, bc.(node.children, Ref(mml)))
end

@inline function pointwise_reduce(x::Array{T,N}, y::Array{T,N})::Array{T,N} where {T<:Real, N}
    x.*y
end

@inline function bc(node::T, mml::Array{Array{S,2},1})::Array{S} where {T, S<:Real}
    mml[node.num]*node.data::Array{S,2}
end

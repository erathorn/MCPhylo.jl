
include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
using BenchmarkTools
Random.seed!(42);


function mpre_order(root::T, traversal::Vector{T})::Vector{T} where T<:Node
    push!(traversal, root)
    if root.nchild != 0
        for ch in root.children
            mpre_order(ch, traversal)
        end
    end # if
    return traversal
end # function pre_order!
function mpre_order(root::T)::Vector{T} where T<:Node
    t::Vector{T} = []
    mpre_order(root, t)
    return t
end # function pre_order!

function get_sisters(root::T) where T<:Node
    mother = root.mother
    t::Vector{T} = []
    for ch in mother.children
        if ch != root
            push!(t, ch)
        end
    end
    return t
end

function GradiantLog2(tree_preorder::Vector{T}, pi_::Number, rates::Array{Float64,1},
                     data, blv) where T<:Node


    Up::Array{Float64,3} = ones(size(data))
    Grad_ll::Array{Float64,1} = zeros(length(tree_preorder))


    @views for node in tree_preorder
        node_ind::Int64 = node.num
        if node.root
            # this is the root
            Up[1,:,node_ind] .= pi_
            Up[2,:,node_ind] .= 1.0-pi_

        else
            for sis in get_sisters(node)
                Up[:,:,node_ind] .*= sis.data
            end
            mother = MCPhylo.get_mother(node)
            @assert any(Up[:,:,mother.num] .!= 1.0)
            Up[:,:,node_ind] .*= Up[:,:,mother.num]

            my_mat::Array{Float64,2} = expo2(pi_, blv[node_ind], 1)
            gradient = sum(Up[:,:,node_ind].*(my_mat'*node.data), dims=1)

            Up[:,:,node_ind] =  node.trprobs*Up[:,:, node_ind]

            gradient ./= sum(Up[:,:,node_ind].*node.data, dims=1)

            Grad_ll[node_ind] = sum(gradient)

            if node.nchild > 0
                @inbounds scaler = sum(Up[:,:,node_ind], dims=1)
                @inbounds Up[:,:,node_ind] ./= scaler
            end # if
        end # if
    end # for

    return Grad_ll[1:end-1]

end # function

function expo2(p, t, r)
    mu::Float64 =  1.0 / (2.0 * p * (1-p))
    return [mu*(p-1)*exp(-mu*t) mu*p*exp(-mu*t);
            mu*(p-1)*-exp(-mu*t) mu*p*-exp(-mu*t)]

end
function expo1(p,t,r)
    return [p*exp(-t*r)/((p - 1)*(-p/(p - 1) + 1)) + 1/(-p/(p - 1) + 1) -p/((p - 1)*(-p/(p - 1) + 1)) - p*exp(-t*r)/((p - 1)*(-p/(p - 1) + 1));            1/(-p/(p - 1) + 1) - exp(-t*r)/(-p/(p - 1) + 1)             -p/((p - 1)*(-p/(p - 1) + 1)) + exp(-t*r)/(-p/(p - 1) + 1)]
end




include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
using BenchmarkTools
Random.seed!(42);
using ForwardDiff

mt, df = make_tree_with_data("local/development.nex"); # load your own nexus file
po = post_order(mt);
for node in po
    node.data = df[:,:,node.num]
    node.scaler = zeros(1,size(node.data, 2))
end
blv = get_branchlength_vector(mt);

rates = ones(1);


pi_ = 0.7
f(x) =  MCPhylo.FelsensteinFunction(po, pi_, rates, df, 3132, x)
#@code_warntype MCPhylo.FelsensteinFunction(po, pi_, rates, df, 3132, blv)
#blv1 = ForwardDiff.Dual.(blv)
#blv1 = [ForwardDiff.Dual(x, ForwardDiff.Partials(tuple(zeros(15)...))) for x in blv]
#dcf=ForwardDiff.GradientConfig(f, blv1)
#ForwardDiff.gradient(f, blv)

using Zygote
gs = Flux.pullback(blv,pi_) do
         MCPhylo.FelsensteinFunction(po, pi_, rates, df, 3132, blv)
       end
g(x,y) =  MCPhylo.FelsensteinFunction(po, x, rates, df, 3132, y)
gs1 = Flux.pullback(g, pi_, blv)

l = gs1[2]
r = gs1[1]
Zg = l(1.0)[2]

using Calculus
cg1 = Calculus.gradient(f, blv)
pi_ = 0.9
function valder_z(g, x,y)
    
    gs1 = Flux.pullback(g, x, y)
    gs1[1], gs1[2](1.0)[2]
end

function valder_c1(g, x,y)
    h(y) = g(x,y)
    cg1 = Calculus.gradient(h,y)
    g(x,y), cg1
end
using BenchmarkTools

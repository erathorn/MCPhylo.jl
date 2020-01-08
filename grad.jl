
include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
using BenchmarkTools
Random.seed!(1234);

mt2, df2 = make_tree_with_data("local/development.nex"); # load your own nexus file
po2 = post_order(mt2);
for node in po2
    node.data = df2[:,:,node.num]
end
blv = get_branchlength_vector(mt2);

rates = ones(1);
f(y) = MCPhylo.FelsensteinFunction(po2, 0.925, rates, df2, 838, y)

x = f(blv)

@benchmark f(blv)

y = f(blv2)

Random.seed!(1234);
mt, df = make_tree_with_data_cu("LangData/Sino-Tibetan.cc.phy.nex"); # load your own nexus file
po = MCPhylo.post_order(mt);
for node in po
    node.data = df[:,:,node.num]
end


blv = MCPhylo.get_branchlength_vector(mt);
rates = ones(1);
g(y) = MCPhylo.FelsensteinFunction(po, 0.996, rates, df, 838, y)
g(blv)


function mf22(blv, f)
        #f(y) = MCPhylo.FelsensteinFunction(po, 0.0996, rates, df, 3132, y)
        f(blv), Calculus.gradient(f, blv, :forward)
end



using Calculus


function mf(blv, f)

    Zygote.refresh()
    v, g = Zygote._pullback(f, blv)
end

function mf2(po, rates, df, blv, f)
        #f(y) = MCPhylo.FelsensteinFunction(po, 0.0996, rates, df, 3132, y)
        f(blv), Calculus.gradient(f, blv, :central)
end



function mf23(po, rates, df, x, f)
    f_x = f(x)
    g = similar(x)
    n = length(x)
    for i = 1:n
        Calculus.@forwardrule x[i] epsilon
        oldx = x[i]
        x[i] = oldx + epsilon
        f_xplusdx = f(x)
        x[i] = oldx
        g[i] = (f_xplusdx - f_x) / epsilon
    end
    return f_x, g
end


function mf3(po, rates, df, y,f)
    #f(y) = MCPhylo.FelsensteinFunction(po, 0.0996, rates, df, 3132, y)
    res = DiffResults.GradientResult(y)
    ReverseDiff.gradient!(res, f, y)
    return DiffResults.value(res), DiffResults.gradient(res)
end


v, g = Zygote.pullback(f, blv)
g(Zygote.sensitivity(v))
Zygote.sensitivity(v)
using Calculus
Calculus.gradient(f, blv)


c = MCPhylo.FelsensteinFunction(po, 0.9605, rates, df, 3132, blvt)
Tracker.grad(blvt)
blvt = param(blv)

r = Flux.Tracker.forward(f,blv)
r[2]([1,1,1])

function val_der(f, y)
    res = DiffResults.GradientResult(y)
    cfg = ForwardDiff.GradientConfig(MyTag, y)
    ForwardDiff.gradient!(res, f, y, cfg, Val{false}())


    return DiffResults.value(res), DiffResults.gradient(res)
end


v, gr = val_der(f, blv)


using Calculus
grc = Calculus.gradient(f, blv)

using Flux
Flux.gradient(f, blv)


#@diff MCPhylo.FelsensteinFunction(po, 0.2705, rates, df, 3132,y1)
Calculus.gradient(f, blv)
using Calculus
v1 = exp(-r*rinc)
m11b = pi_ - (pi_ - 1)*exp(-r*rinc)
m22b = pi_*(exp(-r*rinc) - 1) + 1
m21b = pi_ - pi_* exp(-r*rinc)
m12b = (pi_ - 1)*(exp(-r*rinc) - 1)
pi_ = 0.2
r = 1
rinc =0.5
mm2 = [pi_ - (pi_ - 1)*exp(-r*rinc) pi_ - pi_* exp(-r*rinc); (pi_ - 1)*(exp(-r*rinc) - 1) pi_*(exp(-r*rinc) - 1) + 1]
differentiate("(mm*l_data) .* (mm2*r_data)", :mm2)
differentiate("p - (p-1)*exp(-x)", :x)

ty

exp(MCPhylo.FelsensteinFunction(po, 0.2705, rates, df, 1,blv))


g1 = GradiantLog2(MCPhylo.pre_order(mt), 0.2705, rates, df, 1, blv)


function GradiantLog2(tree_preorder::Vector{Node}, pi_::Number, rates::Array{Float64,1},
                     data, n_c::Int64, blv)#::Array{Any}

    root::Node = tree_preorder[1]
    Up::Array{Float64,3} = ones(size(data))#length(tree_preorder)+1, size(root.data)[1], n_c)
    Grad_ll::Array{Float64,1} = zeros(length(tree_preorder))
    gradient::Vector{Float64} = zeros(n_c)

    for node in tree_preorder
        node_ind::Int64 = node.num
        if node.binary == "1"
            # this is the root
            Up[:,1,node_ind] .= pi_
            Up[:,2,node_ind] .= 1.0-pi_
            #@simd for i in 1:n_c
            #    @inbounds Up[node_ind,1,i] = pi_
            #    @inbounds Up[node_ind,2,i] = 1.0-pi_
            #end # for
        else
            sister::Node = MCPhylo.get_sister(node)
            mother::Node = MCPhylo.get_mother(node)


            sisterdata::Array{Float64} = data[:, :, sister.num]
            #sisterdata ./= sum(sisterdata, dims=1)

            #@inbounds Up[node_ind,:,:] = sisterdata#MCPhylo.pointwise_mat(Up[node_ind,:,:], sisterdata, n_c)
            @inbounds Up[:,:,node_ind] = MCPhylo.pointwise_mat(sisterdata, Up[:,:,mother.num,], n_c)

            nodedata::Array{Float64, 2} = data[:, :, node.num]
            #nodedata ./= sum(nodedata, dims=1)

            for i in 1:n_c
                @inbounds my_mat::Array{Float64,2} = expo2(pi_, node.inc_length, 1)
                @inbounds my_mat2::Array{Float64,2} = expo1(pi_, node.inc_length, 1)

                @inbounds a = nodedata[1, i] * my_mat[1,1] + nodedata[2,i] * my_mat[2,1]
                @inbounds b = nodedata[1, i] * my_mat[1,2] + nodedata[2,i] * my_mat[2,2]

                @inbounds gradient[i] = Up[node_ind,1,i] * a + Up[node_ind,2,i]*b

                @inbounds Up[i,1,node_ind] = Up[i,1, node_ind] * my_mat2[1,1] + Up[i,2,node_ind] * my_mat2[2,1]
                @inbounds Up[i, 2,node_ind] = Up[i,1, node_ind] * my_mat2[1,2] + Up[i,2,node_ind] * my_mat2[2,2]
            end

            d = sum(MCPhylo.pointwise_mat(Up[:,:,node_ind],nodedata, n_c), dims=2)

            @inbounds gradient ./= d[1,:]

            @inbounds Grad_ll[node_ind] = sum(gradient)

            if node.nchild != 0
                @inbounds scaler = sum(Up[:,:,node_ind], dims=1)
                @inbounds Up[:,:,node_ind] ./= scaler
            end # if
        end # if
    end # for

    return Grad_ll

end # function

function exponentiate_binary(pi_::Number, t::Float64, r::Float64)::Array{Float64,2}
    # calculate some auxiliary variables
    ext::Float64 = exp(-t*r)
    ext_::Float64 = 1.0-ext
    pj_::Float64 = 1.0-pi_
    v_::Float64 = ext_*pi_
    w_::Float64 = ext_*pj_

    # return the expontiated matrix
    return [ext+v_ 1-(ext+v_);
            1-(ext+w_) ext+w_]
end

function exponentiate_binary_grad(pi_::Number, t::Float64, r::Float64)::Array{Float64, 2}
        exp2::Float64 = exp(-2*pi_*t*r)*pi_
        return [-exp2 exp2;
                exp2 -exp2]
end

function expo2(p, t, r)
    return [-p*exp(-t*r) p*exp(-t*r);
            (1 - p)*exp(-t*r) (p - 1)*exp(-t*r)]
end
function expo1(p,t,r)
    return [p*exp(-t*r)/((p - 1)*(-p/(p - 1) + 1)) + 1/(-p/(p - 1) + 1) -p/((p - 1)*(-p/(p - 1) + 1)) - p*exp(-t*r)/((p - 1)*(-p/(p - 1) + 1));            1/(-p/(p - 1) + 1) - exp(-t*r)/(-p/(p - 1) + 1)             -p/((p - 1)*(-p/(p - 1) + 1)) + exp(-t*r)/(-p/(p - 1) + 1)]
end

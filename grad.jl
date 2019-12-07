
include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
Random.seed!(1234)
using Calculus

mt, df = make_tree_with_data("local/development.nex") # load your own nexus file
po = MCPhylo.post_order(mt)
blv = MCPhylo.get_branchlength_vector(mt)
rates = ones(1)



f(y) = MCPhylo.FelsensteinFunction(po, 0.2705, rates, deepcopy(df), 3132,y)


#@diff MCPhylo.FelsensteinFunction(po, 0.2705, rates, df, 3132,y1)
Calculus.gradient(f, blv)

differentiate("(p - (p-1)*exp(-x))", :x)
differentiate("p - (p-1)*exp(-x)", :x)


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

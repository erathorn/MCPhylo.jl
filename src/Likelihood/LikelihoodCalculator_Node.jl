

function FelsensteinFunction(
    tree_postorder::Vector{N},
    data::Array{F,4},
    Down::Array{F,4},
    pi_::Array{F},
    trmat::Array{F, 4},
    ptg::Array{F, 4},
)::Tuple{F, Vector{F}} where {N<:GeneralNode, F <: Real}
    ll = FelsensteinFunction(tree_postorder, data, Down, pi_, trmat)
    grv= Vector{F}(undef, size(data,4) - 1)
    pre_order_partial::Array{Float64,4} = similar(data)
    @tturbo for r in axes(pre_order_partial, 3), j in axes(pre_order_partial, 2), i in axes(pre_order_partial,1)
        pre_order_partial[i, j, r, tree_postorder[end].num] = pi_[i]
    end
    
    gradi = Array{Float64, 2}(undef, (size(data, 2), size(data,3)))
    tmp_arr = Array{Float64, 3}(undef, (size(data, 1), size(data, 2), size(data,3)))
    @inbounds for node in reverse(tree_postorder)[2:end]
        mother = get_mother(node)

        @tturbo for r in axes(pre_order_partial, 3), j in axes(pre_order_partial, 2), i in axes(pre_order_partial,1)
            pre_order_partial[i, j, r, node.num] = pre_order_partial[i, j, r,  mother.num]
        end
        
        for c in 1:mother.nchild
            child = mother.children[c]
            if child.num != node.num
                turbo_mul!(pre_order_partial, Down, node.num, child.num)
            end
        end
        mygemmturbo!(tmp_arr, ptg, data, node.num)
        myred!(gradi, pre_order_partial, tmp_arr, node.num)
        mygemmturbo_tr!(tmp_arr, trmat, pre_order_partial, node.num)
        grv[node.num] = comb(tmp_arr, data, gradi, node.num)
        if node.nchild > 0
            bymax!(pre_order_partial, tmp_arr, node.num)
        end
        
    end # for
    ll, grv
end



function FelsensteinFunction(
    tree_postorder::Vector{N},
    data::Array{R,4},
    Down::Array{R,4},
    pi_::Array{R},
    trans_probs::Array{R, 4}
)::R where {N<:GeneralNode{<:Real,<:Integer},R<:Real}

    ll = zero(R)
    c_node::N = tree_postorder[end]


    c1num::Int = 0
    c2num::Int = 0
    pnum::Int = 0

    @views @inbounds for node_ind in eachindex(tree_postorder[1:end-1])

        c_node = tree_postorder[node_ind]
        if c_node.nchild > 0
            
            c1num = c_node.children[1].num
            c2num = c_node.children[2].num
            pnum = c_node.num
            for child in c_node.children
                sum_product_loop!(Down, data, trans_probs, child.num)
                turbo_mul!(data, Down, pnum, child.num)
            end
            ll = by_max!(data, ll, pnum)

        end #if
    end # for
    
    c_node = last(tree_postorder)
    
    for child in c_node.children
        sum_product_loop!(Down, data, trans_probs, child.num) 
        turbo_mul!(data, Down, c_node.num, child.num)
    end

    ll = root_sum(data, pi_, c_node.num, ll)

    ll
end



# function FelsensteinFunction(
#     tree_postorder::Vector{N},
#     pi_::Array{Float64},
#     rates::Vector{Float64},
#     U::Matrix,
#     D::Vector,
#     Uinv::Matrix,
#     mu::Float64,
#     data::Array{Float64,3},
#     substitutionModel::Function,
# ) where {N<:GeneralNode}
    
#     ll = fels_ll(
#         tree_postorder,
#         data_ext,
#         Down,
#         pi_,
#         trans_probs
#     )
#     gr = fels_grad(tree_postorder, data_ext, D, U, Uinv, rates, mu, Down, pi_, blv, trans_probs, blv)
#     return ll, gr
# end


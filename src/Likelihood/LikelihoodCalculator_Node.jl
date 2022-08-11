

function FelsensteinFunction(
    tree_postorder::Vector{N},
    data::Array{F,4},
    Down::Array{F,4},
    pi_::Array{F},
    trmat::Array{F, 4},
    ptg::Array{F, 4},
)::Tuple{F, Vector{F}} where {N<:GeneralNode, F <: Real}
    ll = FelsensteinFunction(tree_postorder, data, Down, pi_, trmat)
    grv = Vector{F}(undef, size(data,4) - 1)
    pre_order_partial::Array{F,4} = similar(data)
    root_num = tree_postorder[end].num
    @tturbo check_empty=false for r in axes(pre_order_partial, 3), j in axes(pre_order_partial, 2), i in axes(pre_order_partial,1)
        pre_order_partial[i, j, r, root_num] = pi_[i]
    end
    
    gradi = Array{F, 2}(undef, (size(data, 2), size(data,3)))
    tmp_arr = Array{F, 3}(undef, (size(data, 1), size(data, 2), size(data,3)))
    @inbounds for node in reverse(tree_postorder)[2:end]
        mother = get_mother(node)
        @tturbo check_empty=false for r in CartesianIndices((axes(pre_order_partial, 1), axes(pre_order_partial, 2), axes(pre_order_partial, 3)))
            pre_order_partial[r, node.num] = pre_order_partial[r,  mother.num]
        end
        c_nums = [c.num for c in mother.children if c.num != node.num]
        turbo_mul!(pre_order_partial, Down, node.num, c_nums)
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

    @views @inbounds for node_ind in eachindex(tree_postorder[1:end-1])

        c_node = tree_postorder[node_ind]
        if c_node.nchild > 0
            comb_sum_product_loop!(Down, data, trans_probs, [c.num for c in c_node.children], c_node.num)
            ll = by_max!(data, ll, c_node.num)
        end #if
    end # for
    
    c_node = last(tree_postorder)
    cnums = [c.num for c in c_node.children]
    comb_sum_product_loop!(Down, data, trans_probs, cnums, c_node.num)
        
    ll = root_sum(data, pi_, c_node.num, ll)
    
    ll
end



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
    root_num = tree_postorder[end].num
    @tturbo for r in axes(pre_order_partial, 3), j in axes(pre_order_partial, 2), i in axes(pre_order_partial,1)
        pre_order_partial[i, j, r, root_num] = pi_[i]
    end
    
    gradi = Array{Float64, 2}(undef, (size(data, 2), size(data,3)))
    tmp_arr = Array{Float64, 3}(undef, (size(data, 1), size(data, 2), size(data,3)))
    @inbounds for node in reverse(tree_postorder)[2:end]
        mother = get_mother(node)
        n_num = node.num
        @tturbo for r in axes(pre_order_partial, 3), j in axes(pre_order_partial, 2), i in axes(pre_order_partial,1)
            pre_order_partial[i, j, r, n_num] = pre_order_partial[i, j, r,  mother.num]
        end
        c_nums = [c.num for c in mother.children if c.num != n_num]
        turbo_mul!(pre_order_partial, Down, n_num, c_nums)
        mygemmturbo!(tmp_arr, ptg, data, n_num)
        myred!(gradi, pre_order_partial, tmp_arr, n_num)
        mygemmturbo_tr!(tmp_arr, trmat, pre_order_partial, n_num)
        grv[n_num] = comb(tmp_arr, data, gradi, n_num)
        if node.nchild > 0
            bymax!(pre_order_partial, tmp_arr, n_num)
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
    pnum::Int = 0

    @views @inbounds for node_ind in eachindex(tree_postorder[1:end-1])

        c_node = tree_postorder[node_ind]
        if c_node.nchild > 0
            pnum = c_node.num
            cnums = [c.num for c in c_node.children]
            comb_sum_product_loop!(Down, data, trans_probs, cnums, pnum)
            ll = by_max!(data, ll, pnum)
            
        end #if
    end # for
    
    c_node = last(tree_postorder)
    cnums = [c.num for c in c_node.children]
    comb_sum_product_loop!(Down, data, trans_probs, cnums, c_node.num)
        
    ll = root_sum(data, pi_, c_node.num, ll)
    
    ll
end

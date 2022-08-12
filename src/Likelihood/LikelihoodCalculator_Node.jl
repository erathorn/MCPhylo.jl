

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
    
    gradi = Array{F, 2}(undef, (size(data, 2), size(data, 3)))
    
    @tturbo for  k in axes(pre_order_partial, 3), j in axes(pre_order_partial, 2), i in axes(pre_order_partial,1)
        pre_order_partial[i, j, k,root_num] = pi_[i]
    end

    for node in reverse(tree_postorder)[2:end]
        mother = get_mother(node)
            
            @tturbo for r in CartesianIndices((axes(pre_order_partial, 1), axes(pre_order_partial, 2), axes(pre_order_partial, 3)))
                pre_order_partial[r, node.num] = pre_order_partial[r, mother.num]
            end
            
            c_nums = [c.num for c in mother.children if c.num != node.num]
            turbo_mul!(pre_order_partial, Down, node.num, c_nums)
            
            @tturbo check_empty=false for n in axes(data, 2), r in axes(data, 3)
                g = zero(eltype(gradi))
                for m ∈ axes(ptg, 1), k ∈ axes(ptg, 2)
                    g += ptg[m, k, r, node.num] * data[k, n, 1, node.num]* pre_order_partial[m, n, r, node.num]
                end
                gradi[n, r] = g
            end
            for r in axes(pre_order_partial, 3)
                pre_order_partial[:, :,  r, node.num] .= transpose(trmat[:, :, r, node.num]) * pre_order_partial[:, :, r, node.num]
            end
            
            res = zero(eltype(grv))
            @tturbo check_empty = false for j in axes(data, 2), r in axes(data, 3)
                col = zero(eltype(grv))
                for k in axes(data, 1)
                    col += data[k, j, r, node.num] * pre_order_partial[k, j, r, node.num]
                end
                res += gradi[j, r] / col
            end
            
            grv[node.num] = res
            
            if node.nchild > 0
                bymax!(pre_order_partial, node.num)    
            end
    end # for
    #@tturbo check_empty=false 
    
    # pre_order_partial[:, :, 1, root_num] .= pi_
    # #
    # #tmp_arr = Array{F, 3}(undef, (size(data, 1), size(data, 2), size(data,3)))
    # @inbounds for node in reverse(tree_postorder)[2:end]
    #     mother = get_mother(node)
    #     #@tturbo check_empty=false 
    #     for r in CartesianIndices((axes(pre_order_partial, 1), axes(pre_order_partial, 2), axes(pre_order_partial, 3)))
    #         pre_order_partial[r, node.num] = Down[r,  mother.num]
    #     end
    #     c_nums = [c.num for c in mother.children if c.num != node.num]
    #     for c in 1:mother.nchild
    #                      child = mother.children[c]
    #                      if child.num != node.num
    #                          pre_order_partial[:, :, node.num] .*= Down[:, :, child.num]
    #                      end
    #                  end
    #     #turbo_mul!(pre_order_partial, Down, node.num, c_nums)
    #     #mygemmturbo!(tmp_arr, ptg, data, node.num)
    #     #myred!(gradi, pre_order_partial, tmp_arr, node.num)
    #     #ptg1 = U * ((D .* (1.0 * mu)) .* diagm(exp.(D .* (1.0*mu*node.inc_length)))) * Uinv
        
    #     gradi = sum(pre_order_partial[:, :,  1, node.num] .* (ptg[:, :, 1, node.num] * data[:, :, 1, node.num]),dims=1)
    #     #gradi = sum(pre_order_partial[:, :,  1, node.num] .* (ptg1 * data[:, :, 1, node.num]),dims=1)
    #     #tr = MCPhylo.calculate_transition(Restriction, 1.0, mu, node.inc_length, U, Uinv, D, pi_)
    #     tr = trmat[:, :, 1, node.num]
    #     pre_order_partial[:, :,  1, node.num] .= transpose(tr) * pre_order_partial[:, :, 1,node.num]
    #     #grv[node.num] = sum(gradi ./ sum(pre_order_partial[:, :, 1, node.num] .* data[:, :,  1,node.num], dims=1))
    #     #mygemmturbo_tr!(tmp_arr, trmat, pre_order_partial, node.num)
    #     grv[node.num] = sum(gradi ./ sum(pre_order_partial[:, :,  1, node.num] .* data[:, :, 1, node.num], dims=1))#comb(tmp_arr, data, gradi, node.num)
    #     if node.nchild > 0
    #         pre_order_partial[:, :,  node.num] ./= maximum(pre_order_partial[:, :,  node.num], dims = 1)
    #         #bymax!(pre_order_partial, tmp_arr, node.num)
    #     end
    

    # end # for
    ll, grv

    # @inbounds @views for node in reverse(tree_postorder)[2:end]
    #     mother = get_mother(node)
    #         pre_order_partial[:, :, node.num] .= pre_order_partial[:, :, mother.num]
    #         @simd for c in 1:mother.nchild
    #             child = mother.children[c]
    #             if child.num != node.num
    #                 pre_order_partial[:, :, node.num] .*= Down[:, :, child.num]
    #             end
    #         end
                
    #         ptg = U * ((D .* (rate * mu)) .* diagm(exp.(D .* (rate*mu*node.inc_length)))) * Uinv
    #         gradi = sum(pre_order_partial[:, :,  node.num] .* (ptg * data[:, :, node.num]),dims=1)
    #         tr = MCPhylo.calculate_transition(substitutionModel, rate, mu, node.inc_length, U, Uinv, D, pi_)
    #         pre_order_partial[:, :,  node.num] .= transpose(tr) * pre_order_partial[:, :, node.num]
    #         grv[node.num] = sum(gradi ./ sum(pre_order_partial[:, :,  node.num] .* data[:, :,  node.num], dims=1))
            
    #         if node.nchild > 0
    #             pre_order_partial[:, :,  node.num] ./= maximum(pre_order_partial[:, :,  node.num], dims = 1)
    #         end
    # end # for
    # grv
end



function FelsensteinFunction(
    tree_postorder::Vector{N},
    data::Array{R,4},
    Down::Array{R,4},
    pi_::Array{R},
    trans_probs::Array{R, 4}
)::R where {N<:GeneralNode{<:Real,<:Integer},R<:Real}

    ll = zero(R)
    
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



function FelsensteinFunction(
    tree_postorder,
    data::Array{F,4},
    Down::Array{F,4},
    pi_::Array{F},
    trmat::Array{F,4},
    ptg::Array{F,4},
)::Tuple{F,Vector{F}} where F<:Real
    ll = FelsensteinFunction(tree_postorder, data, Down, pi_, trmat)
    grv = Vector{F}(undef, size(data, 4) - 1)
    pre_order_partial::Array{F,4} = similar(data)

    gradi = Array{F,2}(undef, (size(data, 2), size(data, 3)))
    tmp = Array{F,2}(undef, (size(data, 1), size(data, 2)))

    @turbo check_empty = false for k in axes(pre_order_partial, 3),
        j in axes(pre_order_partial, 2),
        i in axes(pre_order_partial, 1)

        pre_order_partial[i, j, k, tree_postorder[end].num] = pi_[i]
    end

    for node in reverse(collect(tree_postorder))[2:end]
        mother = get_mother(node)

        @tturbo check_empty = false for r in CartesianIndices((
            axes(pre_order_partial, 1),
            axes(pre_order_partial, 2),
            axes(pre_order_partial, 3),
        ))
            pre_order_partial[r, node.num] = pre_order_partial[r, mother.num]
        end

        @inbounds c_nums = [c.num for c in mother.children if c.num != node.num]
        turbo_mul!(pre_order_partial, Down, node.num, c_nums)

        @tturbo check_empty = false for n in axes(data, 2), r in axes(data, 3)
            g = zero(eltype(gradi))
            for m ∈ axes(ptg, 1), k ∈ axes(ptg, 2)
                g +=
                    ptg[m, k, r, node.num] *
                    data[k, n, r, node.num] *
                    pre_order_partial[m, n, r, node.num]
            end
            gradi[n, r] = g
        end

        @inbounds for r in axes(pre_order_partial, 3)
            mygemmturbo!(
                tmp,
                transpose(trmat[:, :, r, node.num]),
                pre_order_partial[:, :, r, node.num],
            )
            @turbo pre_order_partial[:, :, r, node.num] .= tmp
        end

        res = zero(eltype(grv))
        
        for j in axes(data, 2), r in axes(data, 3)
            col = zero(eltype(grv))
            for k in axes(data, 1)
                col += data[k, j, r, node.num] * pre_order_partial[k, j, r, node.num]
            end
            res += gradi[j, r] / col
        end

        @inbounds grv[node.num] = res

        if node.nchild > 0
            bymax!(pre_order_partial, node.num)
        end
    end # for

    ll, grv

end



function FelsensteinFunction(
    tree_postorder,
    data::Array{R,4},
    Down::Array{R,4},
    pi_::Array{R},
    trans_probs::Array{R,4},
)::R where R<:Real

    ll = zero(R)

    @views @inbounds for node_ind in eachindex(tree_postorder[1:end-1])

        c_node = tree_postorder[node_ind]
        if c_node.nchild > 0
            comb_sum_product_loop!(
                Down,
                data,
                trans_probs,
                [c.num for c in c_node.children],
                c_node.num,
            )
            ll = by_max!(data, ll, c_node.num)
        end #if
    end # for

    c_node = last(tree_postorder)
    @inbounds cnums = [c.num for c in c_node.children]
    comb_sum_product_loop!(Down, data, trans_probs, cnums, c_node.num)

    ll = root_sum(data, pi_, c_node.num, ll)

    ll
end

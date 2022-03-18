
"""
    FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Array{Float64},
                         U::Array{Float64,2}, D::Array{Float64}, Uinv::Array{Float64,2},
                         data::Array{Float64,4}, c_grad::Bool = true) where {N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm. If `c_grad` equals `true` (default) the analytic gradient
regarding the branch lengths of the tree gets computed as well.
"""
function FelsensteinFunction(
    tree_postorder::Vector{N},
    pi_::Array{Float64},
    rates::Float64,
    U::A,
    D::V,
    Uinv::A,
    mu::Float64,
    data::Array{Float64,3},
    substitutionModel::Function,
    c_grad::Bool = true,
) where {N<:GeneralNode,M<:Number,A<:AbstractArray,V<:AbstractVector}
    Nbases, Nsites, Nnodes = size(data)


    grv::Vector{Float64} = Vector{Float64}(undef, Nnodes - 1)
    Down::Array{Float64,3} = similar(data)
    pre_order_partial::Array{Float64,3} = similar(data)
    ll = fels_ll(tree_postorder, data, D, U, Uinv, rates, mu, Down, pi_, substitutionModel)

    if c_grad
        grv = fels_grad(
            tree_postorder,
            data,
            D,
            U,
            Uinv,
            rates,
            mu,
            Down,
            pi_,
            pre_order_partial,
            grv,
            substitutionModel,
        )
    end

    return ll, grv, data, Down, pre_order_partial

end

function fels_grad(
    tree_postorder::Vector{N},
    data::Array{Float64,3},
    D::V,
    U::A,
    Uinv::A,
    rate::Float64,
    mu::Float64,
    Down::Array{Float64,3},
    pi_::Array{Float64},
    pre_order_partial::Array{Float64,3},
    grv::Vector{Float64},
    substitutionModel::Function,
)::Vector{Float64} where {N<:GeneralNode,M<:Number,A<:AbstractArray,V<:AbstractVector}

    root_node = last(tree_postorder)
    pre_order_partial[:, :, root_node.num] .= pi_

    @inbounds @views for node in reverse(tree_postorder)[2:end]
        mother = get_mother(node)
        pre_order_partial[:, :, node.num] .= pre_order_partial[:, :, mother.num]
        for c = 1:mother.nchild
            child = mother.children[c]
            if child.num != node.num
                pre_order_partial[:, :, node.num] .*= Down[:, :, child.num]
            end
        end

        ptg =
            U *
            ((D .* (rate * mu)) .* diagm(exp.(D .* (rate * mu * node.inc_length)))) *
            Uinv
        gradi =
            sum(pre_order_partial[:, :, node.num] .* (ptg * data[:, :, node.num]), dims = 1)
        tr = calculate_transition(
            substitutionModel,
            rate,
            mu,
            node.inc_length,
            U,
            Uinv,
            D,
            pi_,
        )
        pre_order_partial[:, :, node.num] .=
            transpose(tr) * pre_order_partial[:, :, node.num]
        grv[node.num] = sum(
            gradi ./
            sum(pre_order_partial[:, :, node.num] .* data[:, :, node.num], dims = 1),
        )

        if node.nchild > 0
            pre_order_partial[:, :, node.num] ./=
                maximum(pre_order_partial[:, :, node.num], dims = 1)
        end
    end # for
    grv
end

function fels_ll(
    tree_postorder::Vector{N},
    data::Array{Float64,3},
    D::V,
    U::A,
    Uinv::A,
    rate::Float64,
    mu::Float64,
    Down::Array{Float64,3},
    pi_::Array{Float64},
    substitutionModel::Function,
)::Float64 where {N<:GeneralNode,M<:Number,A<:AbstractArray,V<:AbstractVector}



    ll::Float64 = 0.0
    @inbounds @views for node in tree_postorder
        if node.nchild > 0
            data[:, :, node.num] .= 1.0

            for child in node.children
                mul!(
                    Down[:, :, child.num],
                    calculate_transition(
                        substitutionModel,
                        rate,
                        mu,
                        child.inc_length,
                        U,
                        Uinv,
                        D,
                        pi_,
                    ),
                    data[:, :, child.num],
                )

                data[:, :, node.num] .*= Down[:, :, child.num]

            end

            if !node.root
                ll += sum(log.(maximum(data[:, :, node.num], dims = 1)))
                data[:, :, node.num] ./= maximum(data[:, :, node.num], dims = 1)
                #ll += sum(log.(scaler))
            else
                ll += sum(log.(sum(data[:, :, last(tree_postorder).num] .* pi_, dims = 1)))
            end
        end #if
    end # for
    ll
end


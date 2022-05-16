
"""
    FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Array{Float64},
                         U::Array{Float64,2}, D::Array{Float64}, Uinv::Array{Float64,2},
                         data::Array{Float64,4}, c_grad::Bool = true) where {N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm. If `c_grad` equals `true` (default) the analytic gradient
regarding the branch lengths of the tree gets computed as well.
"""
function FelsensteinFunction_wg(
    tree_postorder::Vector{N},
    pi_::Array{Float64},
    rates::Float64,
    U::Matrix,
    D::Vector,
    Uinv::Matrix,
    mu::Float64,
    data::Array{Float64,3},
    substitutionModel::Function,
)::Tuple{Float64, Vector{Float64}, Array{Float64,3},Array{Float64,3},Array{Float64,3}} where N<:GeneralNode
    Nnodes = size(data, 3)


    grv::Vector{Float64} = Vector{Float64}(undef, Nnodes - 1)
    Down::Array{Float64,3} = similar(data)
    pre_order_partial::Array{Float64,3} = similar(data)
    ll = fels_ll(tree_postorder, data, D, U, Uinv, rates, mu, Down, pi_, substitutionModel)

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
    

    return ll, grv, data, Down, pre_order_partial

end


function FelsensteinFunction(
    tree_postorder::Vector{N},
    pi_::Array{Float64},
    rates::Float64,
    U::Matrix,
    D::Vector,
    Uinv::Matrix,
    mu::Float64,
    data::Array{Float64,3},
    substitutionModel::Function,
)::Float64 where N<:GeneralNode
    Down::Array{Float64,3} = similar(data)
    ll = fels_ll(tree_postorder, data, D, U, Uinv, rates, mu, Down, pi_, substitutionModel)
    return ll

end


function fels_grad(
    tree_postorder::Vector{N},
    data::Array{Float64,3},
    D::Vector{R},
    U::Matrix{R},
    Uinv::Matrix{R},
    rate::Float64,
    mu::Float64,
    Down::Array{Float64,3},
    pi_::Array{Float64},
    pre_order_partial::Array{Float64,3},
    grv::Vector{Float64},
    substitutionModel::Function,
)::Vector{Float64} where {N<:GeneralNode,R<:Real}

    @inbounds pre_order_partial[:, :, tree_postorder[end].num] .= pi_
    @inbounds @views for node in reverse(tree_postorder)[2:end]
        mother = get_mother(node)
        pre_order_partial[:, :, node.num] .= pre_order_partial[:, :, mother.num]
        for c in 1:mother.nchild
            child = mother.children[c]
            if child.num != node.num
                pre_order_partial[:, :, node.num] .*= Down[:, :, child.num]
            end
        end
        @fastmath ptg =
            U *
            diagm((D .* (rate * mu)) .* exp.(D .* (rate * mu * node.inc_length))) *
            Uinv
        gradi = redtp(pre_order_partial[:, :, node.num] ,(ptg * data[:, :, node.num]))
        
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
        @fastmath pre_order_partial[:, :, node.num] .=
            transpose(tr) * pre_order_partial[:, :, node.num]
        
        
        grv[node.num] = reddp(gradi, redtp(pre_order_partial[:, :, node.num], data[:, :, node.num]))

        if node.nchild > 0
            bymax!(pre_order_partial[:, :, node.num])
        end
    end # for
    grv
end


function bymax2!(m, ll::F)::F where F<:Real
    maxi = m[1, :]
    ncol = size(m, 2)
    @inbounds for i in axes(m, 1)[2:end]
        for j in 1:ncol
            maxi[j] = maxi[j] < m[i, j] ? m[i, j] : maxi[j]
        end
    end
    @inbounds for i in axes(maxi, 1)
        @fastmath ll += log(maxi[i])
        @fastmath maxi[i] = 1/maxi[i]
    end
    
    @inbounds for i in axes(m, 1)
        for j in axes(m, 2)
            @fastmath m[i, j] *= maxi[j]
        end
    end
    ll
end


function bymax!(m)#::A)::Nothing where A<:AbstractArray{<:Real, 2}
    maxi = m[:, 1]
    ncol = size(m, 2)
    @inbounds for j in 2:ncol
        for i in axes(m, 1)
            maxi[i] = maxi[i] < m[i, j] ? m[i, j] : maxi[i]
        end
    end

    maxi .= 1 ./ maxi
    @inbounds for j in axes(m, 2)
        for i in axes(m, 1)
            @fastmath m[i, j] *= maxi[i]
        end
    end
    nothing
end

function fels_ll(
    tree_postorder::Vector{N},
    data::Array{Float64,3},
    D::Vector{R},
    U::Matrix{R},
    Uinv::Matrix{R},
    rate::Float64,
    mu::Float64,
    Down::Array{Float64,3},
    pi_::Array{Float64},
    substitutionModel::Function,
)::Float64 where {N<:GeneralNode,R<:Real}
    
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
                ll = bymax2!(data[:, :, node.num], ll)
            else
                ll = rootsum(data[:, :, last(tree_postorder).num], pi_, ll)
            end
        end #if
    end # for
    ll
end

function redtp(d1, d2)::Vector
    out = zeros(eltype(d1), axes(d1, 2))
    @inbounds for i in axes(d1, 1)
        for j in axes(d1, 2)
            @fastmath out[j] += d1[i, j] * d2[i, j]
        end
    end
    out
end


function rootsum(data, pi_, ll::F)::F where F
    tmp = zeros(F, axes(data, 2))
    @inbounds for i in axes(data, 1)
        for j in axes(data, 2)
            @fastmath tmp[j] += data[i, j] * pi_[i]
        end
    end
    @inbounds for j in axes(tmp, 1)
        @fastmath ll+= log(tmp[j])
    end
    ll
end

function reddp(d1::A, d2::A) where A
    res = zero(eltype(d1))
    @inbounds for i in axes(d1, 1)
        @fastmath res += d1[i]/d2[i]
    end
    res
end



function fels_ll(
    tree_postorder::Vector{N},
    data::Array{Float64,4},
    D::Vector{R},
    U::Matrix{R},
    Uinv::Matrix{R},
    rates::Vector{Float64},
    mu::Float64,
    pi_::Array{Float64},
    substitutionModel::Function,
    blv::Vector{Float64}
)::Float64 where {N<:GeneralNode{<:Real, <:Integer},R<:Real}
    
    ll::Float64 = 0.0
    c_node::N = tree_postorder[end]
    
    trans_probs = parallel_transition_prob(substitutionModel, U, D, Uinv, rates, mu, blv, pi_)
    tmp_data = deepcopy(data)

    c1num::Int = 0
    c2num::Int = 0
    pnum::Int = 0
    
    @views @inbounds for node_ind in eachindex(tree_postorder[1:end-1])

        c_node = tree_postorder[node_ind]
        if c_node.nchild > 0
                        
            c1num = c_node.children[1].num
            c2num = c_node.children[2].num
            pnum = c_node.num
            
            sum_product_loop!(tmp_data, data, trans_probs, c1num)
            sum_product_loop!(tmp_data, data, trans_probs, c2num)
            turbo_mul!(data, tmp_data, pnum, c1num, c2num)
            ll = bymax22!(data, ll, pnum)
            
        end #if
    end # for

    c_node = last(tree_postorder)
    c1_ind = c_node.children[1].num
    c2_ind = c_node.children[2].num
    c3_ind = c_node.children[3].num
    sum_product_loop!(tmp_data, data, trans_probs, c1_ind)
    sum_product_loop!(tmp_data, data, trans_probs, c2_ind)
    sum_product_loop!(tmp_data, data, trans_probs, c3_ind)
    
    turbo_mul!(data, tmp_data, c_node.num, c1_ind, c2_ind, c3_ind)
    
    
    ll = rootsum2(data[:, :, :, c_node.num], pi_, ll)

    ll
end

function parallel_transition_prob(substitution_model::Function, E::Matrix, D::Vector, EI::Matrix, rates::Vector, mu::Float64, time_vec::Vector, eq_f::Vector)::Array{Float64,4}
    
    out = Array{Float64, 4}(undef, (length(D), length(D),length(rates), length(time_vec),))

     Threads.@threads for b in 1:length(rates)*length(time_vec)
        (r_ind, t_ind) = Base._ind2sub((length(rates), length(time_vec)), b)
        calculate_transition!(substitution_model, out, rates[r_ind],mu, time_vec[t_ind], E, EI, D, eq_f, r_ind, t_ind)
    end
    out
end

function turbo_mul!(data, tmp_data, pnum, c1num, c2num)::Nothing
    @tturbo for s in axes(tmp_data, 1), c in axes(tmp_data, 2), r in axes(tmp_data, 3)
        data[s, c, r, pnum] = tmp_data[s, c, r, c1num] * tmp_data[s, c, r, c2num]
    end
    nothing
end

function turbo_mul!(data, tmp_data, pnum, c1num, c2num, c3num)::Nothing
    @tturbo for s in axes(tmp_data, 1), c in axes(tmp_data, 2), r in axes(tmp_data, 3)
        data[s, c, r, pnum] = tmp_data[s, c, r, c1num] * tmp_data[s, c, r, c2num] * tmp_data[s, c, r, c3num]
    end
    nothing
end


function sum_product_loop!(tmp_data::A, data::A, transprobs::A, num::N)::Nothing where {A, N}
    
    @tturbo for s in axes(data, 1), c in axes(data, 2), r in axes(data, 3)
        tmp = zero(eltype(data))
        for s1 in axes(data, 1)
            tmp += data[s1, c, r, num] * transprobs[s1, s, r, num]
        end
        tmp_data[s, c, r, num] = tmp
    end
    nothing
end



function sum_product_loop!(data, transprobs, pnum, c1num, c2num,c3num, states, columns, rates)::Nothing
    
    @tturbo for s in 1:states, c in 1:columns, r in 1:rates        
        tmp = zero(eltype(data))
        tmp2 = zero(eltype(data))
        tmp3 = zero(eltype(data))
        for s1 in 1:states
            tmp += data[s1, c, r, c1num] * transprobs[s1, s, r, c1num]
            tmp2 += data[s1, c, r, c2num] * transprobs[s1, s, r, c2num]
            tmp3 += data[s1, c, r, c3num] * transprobs[s1, s, r, c3num]
        end
        data[s, c, r, pnum] = tmp * tmp2 * tmp3
    end
    nothing
end


function FelsensteinFunction(
    tree_postorder::Vector{N},
    pi_::Array{Float64},
    rates::Vector{Float64},
    U::Matrix,
    D::Vector,
    Uinv::Matrix,
    mu::Float64,
    data::Array{Float64,3},
    substitutionModel::Function,
)::Float64 where N<:GeneralNode
    blv =get_branchlength_vector(last(tree_postorder))
    data_ext = my_repeat(data, length(rates), [l.num for l in get_leaves(last(tree_postorder))])
    ll = fels_ll(tree_postorder, data_ext, D, U, Uinv, rates, mu, pi_, substitutionModel, blv)
    
    return ll

end


function my_repeat(data::Array{F, N}, nrates::Int, linds::Vector{T})::Array{F, 4} where {F, N, T}
    x, y, z = size(data)
    out = Array{F, 4}(undef, x, y, nrates, z)
    @tturbo for l in eachindex(linds), rind in 1:nrates, xs in axes(out, 1), ys in axes(out, 2)
        out[xs, ys, rind, linds[l]] = data[xs, ys, linds[l]]
    end
    out
end

function bymax22!(m::Array, ll::F, nodenum::Int)::F where F<:Real
    @inbounds maxi = m[1, :, :, nodenum]
    @inbounds for i in axes(m, 1), j in axes(m, 2), k in axes(m,3)
        maxi[j, k] = maxi[j, k] < m[i, j, k, nodenum] ? m[i, j, k, nodenum] : maxi[j,k]     
    end
    @tturbo for i in eachindex(maxi)
        ll += log(maxi[i])
        maxi[i] = 1/maxi[i]
    end
    
    @tturbo for k in axes(m, 3), j in axes(m, 2), i in axes(m, 1)
        m[i, j, k, nodenum] *= maxi[j, k]
    end
    ll
end


function rootsum2(data, pi_, ll::F)::F where F
    
    @tturbo for j in axes(data, 2), k in axes(data, 3)
        tmp = zero(F)
        for i in axes(data, 1)
            tmp += data[i, j, k] * pi_[i]
        end
        ll += log(tmp)
    end

    ll
end


function calculate_transition!(f::typeof(JC), out::A1, rate::R, mu::R, time::R1, U::A, Uinv::A, D::Vector, pi_::Vector, rind, lind)::Nothing where {R1<:Real, R<:Real, A<:AbstractArray{<:Real}, A1<:AbstractArray{<:Real}}
    
    t = rate * time
    if t < MCP_TIME_MIN
        #return_mat = similar(U)
        @inbounds out[:, :, rind, lind] .= 0.0
        @inbounds out[diagind(out[:, :, rind, lind]), rind, lind] .= 1.0
        #return return_mat
    elseif t > MCP_TIME_MAX
        #return_mat = similar(U)
        @inbounds out[:, :, rind, lind] .= 1.0/length(pi_)
        #return return_mat
    else
        t *= mu
        C = similar(U)
         #.=  (U * diagm(exp.(D .* t))) * Uinv
        mygemmturbo!(C, U, diagm(exp.(D .* t)))
        mygemmturbo!(out, C, Uinv, rind, lind)
        #return
    end
    #return_mat
    nothing
end


function mygemmturbo!(C, A, B)
    @tturbo for m ∈ axes(A, 1), n ∈ axes(B, 2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k] * B[k, n]
        end
        C[m, n] = Cmn
    end
end

function mygemmturbo!(C, A, B, rind, lind)
    @tturbo for m ∈ axes(A, 1), n ∈ axes(B, 2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k] * B[k, n]
        end
        C[m, n, rind, lind] = Cmn
    end
end

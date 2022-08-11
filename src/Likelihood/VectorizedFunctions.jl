
function myred!(out::T, d1::R, d2::S, lind::Int)::Nothing where {T, R, S}
    @tturbo check_empty=false for k in eachindex(CartesianIndices(out))
        tmp = zero(eltype(out))
        for i in axes(d1, 1)
            tmp += d1[i, k, lind] * d2[i, k]
        end
        out[k] = tmp
    end
    nothing
end


function mygemmturbo!(C::R, A::T, B::T, lind::Int)::Nothing where {T, R}
    @tturbo check_empty=false for m ∈ axes(A, 1), n ∈ axes(B, 2), r in axes(B, 3)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k, r, lind] * B[k, n, r, lind]
        end
        C[m, n, r] = Cmn
    end
    nothing
end

function mygemmturbo_tr!(C::R, A::T, B::T, lind::Int)::Nothing where {T, R}
    @tturbo check_empty=false for m ∈ axes(A, 1), n ∈ axes(B, 2), r in axes(B, 3)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[k, m, r, lind] * B[k, n, r, lind]
        end
        C[m, n, r] = Cmn
    end
    nothing
end




function comb(po::C, data::Array{B, 4}, gradi::A, lind::Int)::B where {A, B, C}
    res = zero(eltype(data))
    @tturbo check_empty=false for r in eachindex(CartesianIndices(gradi))#k in axes(data, 3), j in axes(data, 2)
        tmp = zero(eltype(data))
        for i in axes(data, 1)
            tmp += data[i, r, lind] * po[i, r]
        end
        res += gradi[r] / tmp
    end
    res
end


function my_repeat(
    data::Array{F,N},
    nrates::Int,
    linds::Vector{T},
    nlinds::Vector{T},
)::Array{F,4} where {F,N,T}
    x, y, z = size(data)
    out = Array{F,4}(undef, x, y, nrates, z)
    @tturbo check_empty=false for l in eachindex(linds),
        rind in axes(out, 3),
        xs in axes(out, 1),
        ys in axes(out, 2)

        out[xs, ys, rind, linds[l]] = data[xs, ys, linds[l]]
    end
    @tturbo check_empty=false for l in eachindex(nlinds),
        rind in axes(out, 3),
        xs in axes(out, 1),
        ys in axes(out, 2)

        out[xs, ys, rind, nlinds[l]] = one(eltype(out))
    end
    out
end

function by_max!(m::Array, ll::F, nodenum::Int)::F where {F<:Real}
    @inbounds maxi = m[1, :, :, nodenum]
    @inbounds for k in axes(m, 3), j in axes(m, 2), i in axes(m, 1)[2:end]
        maxi[j, k] = maxi[j, k] < m[i, j, k, nodenum] ? m[i, j, k, nodenum] : maxi[j, k]
    end
    @tturbo check_empty=false for i in eachindex(maxi)
        ll += log(maxi[i])
        maxi[i] = 1 / maxi[i]
    end

    @tturbo check_empty=false for i in axes(m, 1), k in eachindex(CartesianIndices(maxi))
        m[i, k, nodenum] *= maxi[k]
    end
    ll
end


function bymax!(out::A, m::B, nodenum::Int)::Nothing where {A, B}
    @inbounds maxi = m[:, 1, :]
    @inbounds for r in axes(m, 3), j in axes(m, 2)[2:end], i in axes(m, 1)
        maxi[i, r] = maxi[i, r] < m[i, j, r] ?  m[i, j, r] : maxi[i, r]
    end
    
    @tturbo check_empty=false for r in axes(m, 3), j in axes(m, 2), i in axes(m, 1)      
        out[i, j, r, nodenum] = m[i, j, r] / maxi[i, r]
    end
    nothing
end


function root_sum(data::T, pi_::S, lind::Int, ll::F)::F where {F, T, S}

    @tturbo check_empty=false for k in CartesianIndices((axes(data, 2), axes(data, 3)))
        tmp = zero(F)
        for i in axes(data, 1)
            tmp += data[i, k, lind] * pi_[i]
        end
        ll += log(tmp)
    end

    ll
end


@inline function diagonalizer!(out::B, D::A, blv::A, mu::C, rates::A)::Nothing where {A, B, C}
    @tturbo check_empty=false for i in eachindex(blv), r in eachindex(rates), d in eachindex(D)
        out[d, d, r, i] = exp(mu * blv[i] * D[d] * rates[r])
    end
    nothing
end

@inline function diagonalizer_ptg!(out::B, D::A, blv::A, mu::C, rates::A)::Nothing where {A, B, C}
    @tturbo check_empty=false for i in eachindex(blv), r in eachindex(rates), d in eachindex(D)
        out[d, d, r, i] = D[d] * rates[r] * mu * exp(mu * blv[i] * D[d] * rates[r])
    end
    nothing
end


function R_gemmturbo_large!(C::T, A::S, B::T)::Nothing where {T, S}
    @tturbo check_empty=false for m ∈ axes(A, 1), l ∈ CartesianIndices((axes(B, 2), axes(B, 3), axes(B, 4)))
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k] * B[k, l]
        end
        C[m, l] = Cmn
    end
    nothing
end


function L_gemmturbo_large!(C::T, A::T, B::S)::Nothing where {T, S}
    @tturbo check_empty=false for r ∈ axes(A, 3), l ∈ axes(A, 4),  m ∈ axes(A, 1), n ∈ axes(B, 2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k, r, l] * B[k, n]
        end
        C[m, n, r, l] = Cmn
    end
    nothing
end


function parallel_transition_prob(
    U::Matrix,
    D::Vector,
    Uinv::Matrix,
    rates::Vector,
    mu::Float64,
    blv::Vector{R},
)::Array{R,4} where {R<:Real}
    out = zeros(R, (length(D), length(D), length(rates), length(blv)))
    out2 = Array{R,4}(undef, (length(D), length(D), length(rates), length(blv)))
    diagonalizer!(out, D, blv, mu, rates)
    R_gemmturbo_large!(out2, U, out)
    L_gemmturbo_large!(out, out2, Uinv)
    out
end

function turbo_mul!(data::A, tmp_data::A, pnum::D, nums::Vector{D})::Nothing where {A,D}
    @tturbo check_empty=false for num_ind in eachindex(nums), r in CartesianIndices((axes(tmp_data, 1), axes(tmp_data, 2), axes(tmp_data,3)))
        data[r, pnum] *= tmp_data[r, nums[num_ind]]
    end
    nothing
end


function sum_product_loop!(tmp_data::A, data::A, transprobs::A, num::N)::Nothing where {A,N}

    @tturbo check_empty=false for s in axes(data, 1), c in axes(data, 2), r in axes(data, 3)
        tmp = zero(eltype(data))
        for s1 in axes(data, 1)
            tmp += data[s1, c, r, num] * transprobs[s, s1, r, num]
        end
        tmp_data[s, c, r, num] = tmp
    end
    nothing
end


function turbo_dot(x::Vector{T}, y::Vector{T})::T where {T}
    res = zero(T)
    @turbo check_empty=false for i in eachindex(x, y)
        res += x[i] * y[i]
    end
    res
end



function comb_sum_product_loop!(
    tmp_data::A,
    data::A,
    transprobs::A,
    nums::Vector{N},
    pnum::N,
)::Nothing where {A,N}

    @tturbo check_empty=false inline=true for num_ind in eachindex(nums),
        s in axes(data, 1),
        c in axes(data, 2),
        r in axes(data, 3)

        tmp = zero(eltype(data))
        num = nums[num_ind]
        mul_red = data[s, c, r, pnum]
        for s1 in axes(data, 1)
            tmp += data[s1, c, r, num] * transprobs[s, s1, r, num]
        end
        tmp_data[s, c, r, num] = tmp
        data[s, c, r, pnum] = mul_red * tmp
    end
    nothing

end

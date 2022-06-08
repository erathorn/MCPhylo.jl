
function myred!(out, d1, d2, lind)
    @turbo for j in axes(d1, 2), r in axes(d1, 3)
        tmp = zero(eltype(out))
        for i in axes(d1, 1)
            tmp += d1[i, j, r, lind] * d2[i, j, r]
        end
        out[j, r] = tmp
    end
end


function mygemmturbo!(C, A, B, lind)
    @turbo for m ∈ axes(A, 1), n ∈ axes(B, 2), r in axes(B, 3)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k, r, lind] * B[k, n, r, lind]
        end
        C[m, n, r] = Cmn
    end
end

function mygemmturbo_tr!(C, A, B, lind)
    @turbo for m ∈ axes(A, 1), n ∈ axes(B, 2), r in axes(B, 3)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[k, m, r, lind] * B[k, n, r, lind]
        end
        C[m, n, r] = Cmn
    end
end




function comb(po, data, gradi, lind)
    res = zero(eltype(data))
    @turbo for k in axes(data,3), j in axes(data,2)
        tmp = zero(eltype(data))
        for i in axes(data, 1)
            tmp += data[i, j, k, lind] * po[i, j, k]
        end
        res += gradi[j, k] / tmp
    end
    res
end


function my_repeat(
    data::Array{F,N},
    nrates::Int,
    linds::Vector{T},
    nlinds::Vector{T}
)::Array{F,4} where {F,N,T}
    x, y, z = size(data)
    out = Array{F,4}(undef, x, y, nrates, z)
    @tturbo for l in eachindex(linds),
        rind = 1:nrates,
        xs in axes(out, 1),
        ys in axes(out, 2)
        out[xs, ys, rind, linds[l]] = data[xs, ys, linds[l]]
    end
    @tturbo for l in eachindex(nlinds),
        rind = 1:nrates,
        xs in axes(out, 1),
        ys in axes(out, 2)
        out[xs, ys, rind, nlinds[l]] = one(eltype(out))
    end
    out
end

function by_max!(m::Array, ll::F, nodenum::Int)::F where {F<:Real}
    @inbounds maxi = m[1, :, :, nodenum]
    @inbounds for k in axes(m, 3),j in axes(m, 2), i in axes(m, 1)[2:end]
        maxi[j, k] = maxi[j, k] < m[i, j, k, nodenum] ? m[i, j, k, nodenum] : maxi[j, k]
    end
    @turbo for i in eachindex(maxi)
        ll += log(maxi[i])
        maxi[i] = 1 / maxi[i]
    end

    @turbo for k in axes(m, 3), j in axes(m, 2), i in axes(m, 1)
        m[i, j, k, nodenum] *= maxi[j, k]
    end
    ll
end


function bymax!(out, m, nodenum)
    @inbounds maxi = m[:, 1, :]
    @inbounds for r in axes(m,3),j in axes(m,2)[2:end],  i in axes(m, 1)
        maxi[i, r] = maxi[i, r] < m[i, j, r] ? m[i, j, r] : maxi[i, r]
    end

    @turbo for r in axes(m, 3) , i in axes(m,1)
        tmp = 1 / maxi[i, r]
        for j in axes(m, 2)
            out[i, j, r, nodenum] = m[i, j, r]*tmp
        end
    end
    nothing
end


function root_sum(data, pi_, lind, ll::F)::F where {F}

    @turbo for k in axes(data, 3), j in axes(data, 2)
        tmp = zero(F)
        for i in axes(data, 1)
            tmp += data[i, j, k, lind] * pi_[i]
        end
        ll += log(tmp)
    end

    ll
end


function diagonalizer!(out, D, blv, mu, rates)::Nothing
    @turbo for i in eachindex(blv), r in eachindex(rates), d in eachindex(D)
        out[d, d, r, i] = exp(mu * blv[i] * D[d] * rates[r])
    end
    nothing
end

function diagonalizer_ptg!(out, D, blv, mu, rates)::Nothing
    @turbo for i in eachindex(blv), r in eachindex(rates), d in eachindex(D)
        out[d, d, r, i] = D[d] * rates[r] * mu * exp(mu * blv[i] * D[d] * rates[r])
    end
    nothing
end


function R_gemmturbo_large!(C, A, B)
    @tturbo for m ∈ axes(A, 1), l ∈ axes(B, 4), r ∈ axes(B, 3), n ∈ axes(B, 2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k] * B[k, n, r, l]
        end
        C[m, n, r, l] = Cmn
    end
end


function L_gemmturbo_large!(C, A, B)
    @tturbo for l ∈ axes(A, 4), r ∈ axes(A, 3), m ∈ axes(A, 1), n ∈ axes(B, 2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k, r, l] * B[k, n]
        end
        C[m, n, r, l] = Cmn
    end
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

function turbo_mul!(data::A, tmp_data::A, pnum::D, c1num::D, c2num::D)::Nothing where {A, D}
    @turbo for r in axes(tmp_data, 3), c in axes(tmp_data, 2), s in axes(tmp_data, 1)
        data[s, c, r, pnum] = tmp_data[s, c, r, c1num] * tmp_data[s, c, r, c2num]
    end
    nothing
end


function turbo_mul!(data::A, tmp_data::A, pnum::D, c1num::D)::Nothing where {A, D}
    @turbo for r in axes(tmp_data, 3), c in axes(tmp_data, 2), s in axes(tmp_data, 1)
        data[s, c, r, pnum] *= tmp_data[s, c, r, c1num]
    end
    nothing
end


function turbo_mul!(data::A, tmp_data::A, pnum::D, c1num::D, c2num::D, c3num::D)::Nothing where {A, D}
    @turbo for s in axes(tmp_data, 1), c in axes(tmp_data, 2), r in axes(tmp_data, 3)
        data[s, c, r, pnum] =
            tmp_data[s, c, r, c1num] * tmp_data[s, c, r, c2num] * tmp_data[s, c, r, c3num]
    end
    nothing
end


function sum_product_loop!(tmp_data::A, data::A, transprobs::A, num::N)::Nothing where {A,N}

    @tturbo for s in axes(data, 1), c in axes(data, 2), r in axes(data, 3)
        tmp = zero(eltype(data))
        for s1 in axes(data, 1)
            tmp += data[s1, c, r, num] * transprobs[s, s1, r, num]
        end
        tmp_data[s, c, r, num] = tmp
    end
    nothing
end


function turbo_dot(x::Vector{T}, y::Vector{T})::T where T
    res = zero(T)
    @turbo for i in eachindex(x)
        res += x[i] * y[i]
    end
    res
end

function refraction!(
    s::Tree_HMC_State{T},
    epsilon::Float64,
    logfgrad::Function,
    logf::Function,
    delta::Float64,
)::Tuple{Int,Int} where {T<:GeneralNode}

    blenvec = get_branchlength_vector(s.x)

    @. s.r += (epsilon * 0.5) * s.g

    tmpB = @. blenvec + (epsilon * s.r)

    tmpB, nni, att_nni = ref_NNI!(s, tmpB, epsilon, blenvec, delta, logf)

    blenvec = molifier.(tmpB, delta)

    set_branchlength_vector!(s.x, blenvec)

    logf, grad = logfgrad(s.x)

    fac = scale_fac.(blenvec, delta)
    grad .*= fac
    @. s.r += (epsilon * 0.5) * grad
    s.g[:] .= grad[:]
    s.lf = logf

    return nni, att_nni
end



function ref_NNI!(
    s::Tree_HMC_State{T},
    tmpB::V,
    epsilon::Float64,
    blv::V,
    delta::Float64,
    logf::Function,
)::Tuple{V,Int,Int} where {V<:AbstractVector{Float64},T}

    intext = internal_external(s.x)
    t = 0.0
    nni = zero(Int)
    pm = epsilon > 0 ? 1 : -1
    att_nni = zero(Int)
    while minimum(tmpB) <= 0.0

        timelist = tmpB ./ abs.(s.r)
        ref_index = argmin(timelist)

        temp = epsilon - t + timelist[ref_index]
        blv = @. blv + temp * s.r
        blv = abs.(blv)
        s.r[ref_index] *= -1.0

        if intext[ref_index] == 1

            set_branchlength_vector!(s.x, molifier.(blv, delta))
            U_before_nni = logf(s.x)

            v_copy = deepcopy(s.x)
            tmp_NNI_made = NNI!(v_copy, ref_index)

            if tmp_NNI_made != 0
                att_nni += 1
                U_after_nni = logf(v_copy)
                #U_before_nni = fetch(temp)
                delta_U = 2.0 * (U_before_nni - U_after_nni)
                my_v = s.r[ref_index]^2

                if my_v > delta_U
                    nni += tmp_NNI_made
                    s.r[ref_index] = sqrt(my_v - delta_U)
                    s.x = v_copy
                end # if my_v
            #else
            #    _ = fetch(temp) # fetch to free handle
            end #if NNI

        end #non leave

        t = epsilon + timelist[ref_index]
        tmpB = @. blv + (epsilon - t) * s.r

    end #while

    tmpB, nni, att_nni
end

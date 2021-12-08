
function refraction!(
    s::Tree_HMC_State,
    epsilon::Float64,
    logfgrad::Function,
    delta::Float64,
)

    blenvec = get_branchlength_vector(s.x)
    fac = scale_fac.(blenvec, delta)

    @. s.r += (epsilon * 0.5) * s.g * fac

    tmpB = @. blenvec + (epsilon * s.r)

    nni = 0

    if minimum(tmpB) <= 0
        tmpB, nni = ref_NNI!(s, tmpB, abs(epsilon), blenvec, delta, logfgrad)
    end

    blenvec = molifier.(tmpB, delta)

    set_branchlength_vector!(s.x, blenvec)

    logf, grad = logfgrad(s.x)

    fac = scale_fac.(blenvec, delta)

    @. s.r += (epsilon * 0.5) * grad * fac
    s.g[:] .= grad[:]
    s.lf = logf

    return nni
end



function ref_NNI!(
    s::Tree_HMC_State,
    tmpB::Vector{Float64},
    epsilon::Float64,
    blv::Vector{Float64},
    delta::Float64,
    logfgrad::Function,
)

    intext = internal_external(s.x)
    t = 0.0
    nni = 0
    pm = epsilon > 0 ? 1 : -1

    while minimum(tmpB) <= 0.0

        timelist = tmpB ./ abs.(s.r)
        ref_index = argmin(timelist)

        temp = abs(epsilon) - t + timelist[ref_index]
        blv = @. blv + temp * s.r

        s.r[ref_index] *= -1.0

        if intext[ref_index] == 1

            set_branchlength_vector!(s.x, molifier.(blv, delta))

            U_before_nni, _ = logfgrad(s.x) # still with molified branch length

            v_copy = deepcopy(s.x)

            tmp_NNI_made = NNI!(v_copy, ref_index)

            if tmp_NNI_made != 0

                U_after_nni, _ = logfgrad(v_copy)

                delta_U = 2.0 * (U_before_nni - U_after_nni)
                my_v = s.r[ref_index]^2

                if my_v > delta_U
                    nni += tmp_NNI_made
                    s.r[ref_index] = sqrt(my_v - delta_U)
                    s.x = v_copy
                end # if my_v
            end #if NNI
        end #non leave

        t = abs(epsilon) + timelist[ref_index]
        tmpB = @. blv + (pm * (abs(epsilon) - t)) * s.r

    end #while

    tmpB, nni
end

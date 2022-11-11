

function nuts_sub!(
    v::Sampler{NUTSTune{riemann_nuts,F,G},T},
    epsilon::Float64,
    logfgrad::Function,
) where {T<:AbstractArray{<:Real},F,G}
    x = v.value

    r = randn(size(x))

    log_sum_weight = 0.0

    logf, grad = logfgrad(x)


    z_plus = Array_HMC_State(x[:], r[:], grad[:], logf)
    z_minus = Array_HMC_State(x[:], r[:], grad[:], logf)
    z_propose = Array_HMC_State(x[:], r[:], grad[:], logf)
    z_sample = Array_HMC_State(x[:], r[:], grad[:], logf)

    p_fwd_fwd = r[:]
    p_fwd_bck = r[:]
    p_bck_fwd = r[:]
    p_bck_bck = r[:]
    p_sharp_fwd_fwd = r[:]
    p_sharp_fwd_bck = r[:]
    p_sharp_bck_fwd = r[:]
    p_sharp_bck_bck = r[:]

    # calculate hamiltonian of current state
    H0 = hamiltonian(z_sample)

    rho = r

    tnni = 0
    depth = 0


    acc_p_r = 0

    meta = NUTSMeta()

    while depth < v.tune.tree_depth
        meta.alpha = 0
        meta.nalpha = 0


        log_sum_weight_subtree = [-Inf]
        valid_subtree = false

        pm = 2 * (rand() > 0.5) - 1

        rho_bwd = zeros(size(x))
        rho_fwd = zeros(size(x))

        if pm == -1
            z_worker = Ref(z_minus)
            rho_fwd .= rho
            p_fwd_bck .= p_bck_bck
            p_sharp_fwd_bck .= p_sharp_bck_bck

            valid_subtree = buildtree(
                z_worker,
                z_propose,
                pm,
                depth,
                epsilon,
                logfgrad,
                H0,
                rho_bwd,
                log_sum_weight_subtree,
                p_sharp_bck_fwd,
                p_sharp_bck_bck,
                p_bck_fwd,
                p_bck_bck,
                meta,
            )

        else

            z_worker = Ref(z_plus)
            rho_bwd .= rho
            p_bck_fwd .= p_fwd_fwd
            p_sharp_bck_fwd .= p_sharp_fwd_fwd

            valid_subtree = buildtree(
                z_worker,
                z_propose,
                pm,
                depth,
                epsilon,
                logfgrad,
                H0,
                rho_fwd,
                log_sum_weight_subtree,
                p_sharp_fwd_bck,
                p_sharp_fwd_fwd,
                p_fwd_bck,
                p_fwd_fwd,
                meta,
            )

        end#if pm
        tnni += meta.nni
        # the subtree is not valid! Further computations are not necessary
        if !valid_subtree
            break
        end
        depth += 1
        lsw_t = log_sum_weight_subtree[1]

        if lsw_t > log_sum_weight
            transfer!(z_sample, z_propose)
            acc_p_r += 1
        else
            accp = exp(lsw_t - log_sum_weight)
            if rand() < accp
                transfer!(z_sample, z_propose)
                acc_p_r += 1
            end
        end
        log_sum_weight = logaddexp(lsw_t, log_sum_weight)

        rho = rho_bwd .+ rho_fwd

        pers = nouturn(rho, p_sharp_bck_bck, p_sharp_fwd_fwd)

        rho_extend = rho_bwd .+ p_fwd_bck

        pers &= nouturn(rho_extend, p_sharp_bck_bck, p_sharp_fwd_bck)

        rho_extend = rho_fwd .+ p_bck_fwd

        pers &= nouturn(rho_extend, p_sharp_bck_fwd, p_sharp_fwd_fwd)

        if !pers
            break
        end

    end


    v.tune.stepsizeadapter.metro_acc_prob = meta.alpha / meta.nalpha

    v.value[:] .= z_sample.x[:]

    push!(v.tune.acc_p_r, acc_p_r)
    v
end


function buildtree(
    s_worker::Ref{T},
    s_prob::T,
    pm::Int64,
    depth::Integer,
    epsilon::Float64,
    logfgrad::Function,
    H0::R,
    rho::Vector{Float64},
    log_sum_weight::Vector{Float64},
    p_sharp_beg::Vector{Float64},
    p_sharp_end::Vector{Float64},
    p_beg::Vector{Float64},
    p_end::Vector{Float64},
    meta::NUTSMeta,
) where {T<:Array_HMC_State,R<:Real}

    if depth == 0
        leapfrog!(s_worker[], pm * epsilon, logfgrad)
        transfer!(s_prob, s_worker[])


        h = hamiltonian(s_prob)
        h = isnan(h) ? -Inf : h


        divergent = (h - H0) < -1000

        stat = (h - H0) > 0 ? 1 : exp(h - H0)

        meta.alpha += stat
        meta.nalpha += 1

        log_sum_weight .= logaddexp(log_sum_weight[1], h - H0)

        rho .+= s_prob.r
        p_sharp_beg .= s_prob.r
        p_sharp_end .= s_prob.r
        p_beg .= s_prob.r
        p_end .= s_prob.r

        return !divergent
    end

    rho_init = zeros(size(rho))
    log_sum_weight_init = [-Inf]
    p_init_end = zeros(size(rho))
    p_sharp_init_end = zeros(size(rho))

    valid_init = buildtree(
        s_worker,
        s_prob,
        pm,
        depth - 1,
        epsilon,
        logfgrad,
        H0,
        rho_init,
        log_sum_weight_init,
        p_sharp_beg,
        p_sharp_init_end,
        p_beg,
        p_init_end,
        meta,
    )

    if !valid_init
        return false
    end

    s_final = transfer(s_worker[])
    rho_final = zeros(size(rho))
    log_sum_weight_final = [-Inf]
    p_final_beg = zeros(size(rho))
    p_sharp_final_beg = zeros(size(rho))

    valid_final = buildtree(
        s_worker,
        s_final,
        pm,
        depth - 1,
        epsilon,
        logfgrad,
        H0,
        rho_final,
        log_sum_weight_final,
        p_sharp_final_beg,
        p_sharp_end,
        p_final_beg,
        p_end,
        meta,
    )

    if !valid_final
        return false
    end

    log_sum_weight_subtree = logaddexp(log_sum_weight_init[1], log_sum_weight_final[1])
    log_sum_weight[1] = logaddexp(log_sum_weight[1], log_sum_weight_subtree)

    if log_sum_weight_final[1] > log_sum_weight_subtree
        transfer!(s_prob, s_final)
    else
        accp = exp(log_sum_weight_final[1] - log_sum_weight_subtree)
        if rand() < accp
            transfer!(s_prob, s_final)
        end
    end

    rho_subtree = rho_init + rho_final

    rho .+= rho_subtree

    pers_crit = nouturn(rho_subtree, p_sharp_end, p_sharp_beg)

    rho_subtree = rho_init + p_final_beg

    pers_crit &= nouturn(rho_subtree, p_sharp_beg, p_sharp_final_beg)

    rho_subtree = rho_final + p_init_end

    pers_crit &= nouturn(rho_subtree, p_sharp_init_end, p_sharp_end)

    return pers_crit

end


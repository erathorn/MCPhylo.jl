
#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################

mutable struct PNUTS_Rie_Tune <: SamplerTune
    logfgrad::Union{Function,Missing}
    stepsizeadapter::PNUTSstepadapter
    adapt::Bool
    epsilon::Float64
    delta::Float64
    moves::Vector{Int}
    tree_depth::Int    
    tree_depth_trace::Vector{Int}
    acc_p_r::Vector{Int}

    PNUTS_Rie_Tune() = new()

    function PNUTS_Rie_Tune(
        x::Vector{T},
        epsilon::Float64,
        logfgrad::Union{Function,Missing};
        target::Real = 0.75,
        tree_depth::Int = 10,
        targetNNI::Float64 = 5.0,
        delta::Float64 = 0.003
    ) where {T<:GeneralNode}
        
        new(
            logfgrad,
            PNUTSstepadapter(0,0,0,PNUTS_StepParams(0.5,target,0.05,0.75,10,targetNNI)),
            false,
            epsilon,
            delta,
            Int[],
            tree_depth,
            Int[],
            Int[]
        )
    end
end

PNUTS_Rie_Tune(
    x::Vector{T},
    logfgrad::Function,
    ::NullFunction,
    delta::Float64 = 0.003,
    target::Real = 0.75;
    args...,
) where {T<:GeneralNode} =
    PNUTS_Rie_Tune(x, nutsepsilon(x[1], logfgrad, delta, target), logfgrad; args...)

PNUTS_Rie_Tune(
    x::Vector{T},
    logfgrad::Function,
    delta::Float64,
    target::Real;
    args...,
) where {T<:GeneralNode} =
    PNUTS_Rie_Tune(x, nutsepsilon(x[1], logfgrad, delta, target), logfgrad; args...)

PNUTS_Rie_Tune(x::Vector; epsilon::Real, args...) = PNUTS_Rie_Tune(x, epsilon, missing, args...)

const PNUTS_Rie_Variate = SamplerVariate{PNUTS_Rie_Tune}


#################### Sampler Constructor ####################


"""
    PNUTS(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for PNUTS sampling. The Parameter is assumed to be
a tree.

Returns a `Sampler{PNUTSTune}` type object.

* params: stochastic node to be updated with the sampler.

* args...: additional keyword arguments to be passed to the PNUTSVariate constructor.
"""
function PNUTS_Rie(params::ElementOrVector{Symbol}; args...)
    samplerfx = function (model::Model, block::Integer)
        block = SamplingBlock(model, block, true)

        f = let block = block
            (x, sz, ll, gr) -> mlogpdfgrad!(block, x, sz, ll, gr)
        end
        v = SamplerVariate(block, f, NullFunction(); args...)
        
        sample!(v::PNUTS_Rie_Variate, f, adapt = model.iter <= model.burnin)
        
        relist(block, v)
    end
    Sampler(params, samplerfx, PNUTS_Rie_Tune())
end


#################### Sampling Functions ####################

sample!(v::PNUTS_Rie_Variate; args...) = sample!(v, v.tune.logfgrad; args...)

function sample!(v::PNUTS_Rie_Variate, logfgrad::Function; adapt::Bool = false)
    tune = v.tune
    adapter = tune.stepsizeadapter
    const_params = tune.stepsizeadapter.params
    
    setadapt!(v, adapt)
    
    if tune.adapt
        adapter.m += 1
        
        nuts_sub!(v, tune.epsilon, logfgrad)

        adaptstat = adapter.metro_acc_prob > 1 ? 1 : adapter.metro_acc_prob
        
        adaptstat = const_params.δ - adaptstat
        #@show adaptstat, const_params.δ
        HT2 = const_params.τ - adapter.avg_nni
        adaptstat -= HT2#0.4*abs_adapter(HT2)
        #HT2 = tune.targetNNI - tune.nniprime / tune.nalpha
        adaptstat /= 2
        #@show adapter.avg_nni,adapter.metro_acc_prob, adaptstat, HT2#, adaptstat2
        η = 1.0/(adapter.m + const_params.t0)
        
        adapter.s_bar = (1.0 - η) * adapter.s_bar + η * adaptstat
        x = const_params.μ - adapter.s_bar * sqrt(adapter.m) / const_params.γ
        #@show exp(x), adapter.s_bar, η
        x_η = adapter.m^-const_params.κ
        adapter.x_bar = (1.0 - x_η) * adapter.x_bar + x_η * x
        tune.epsilon = exp(x)

    else
        if (adapter.m > 0)
            tune.epsilon = exp(adapter.x_bar)
        end

        nuts_sub!(v, tune.epsilon, logfgrad)
    end
    v
end


function setadapt!(v::PNUTS_Rie_Variate, adapt::Bool)
    tune = v.tune
    if adapt && !tune.adapt
        tune.stepsizeadapter.m = 0
        tune.stepsizeadapter.params = update_step(tune.stepsizeadapter.params, log(10.0 * tune.epsilon))
    end
    tune.adapt = adapt
    v
end



function nuts_sub!(v::PNUTS_Rie_Variate, epsilon::Float64, logfgrad::Function)
    x = deepcopy(v.value[1])
    nl = size(x)[1] - 1
    delta = v.tune.delta
    r = randn(nl)
    #g = zeros(nl)
    #epsilon = 0.15
    #@show epsilon
    
    log_sum_weight = 0.0

    blv = get_branchlength_vector(x)
    set_branchlength_vector!(x, molifier.(blv, delta))
    logf, grad = logfgrad(x, nl, true, true)
    

    z_plus = Tree_HMC_State(deepcopy(x), r, grad, logf)
    z_minus = Tree_HMC_State(deepcopy(x), r, grad, logf)
    z_propose = Tree_HMC_State(deepcopy(x), r, grad, logf)
    z_sample = Tree_HMC_State(deepcopy(x), r, grad, logf)
    z_worker = Tree_HMC_State(deepcopy(x), r, grad, logf)

    p_fwd_fwd = r
    p_fwd_bck = r
    p_bck_fwd = r
    p_bck_bck = r
    p_sharp_fwd_fwd = grad
    p_sharp_fwd_bck = grad
    p_sharp_bck_fwd = grad
    p_sharp_bck_bck = grad

    # calculate hamiltonian of current state
    H0 = hamiltonian(z_sample)
    
    rho = r
    
    nni = 0
    tnni = 0
    depth = 0
        
    
    acc_p_r = 0

    meta = PNUTSMeta()

    while depth < v.tune.tree_depth
        #reset!(meta)
        meta.alpha = 0
        meta.nalpha = 0
        meta.accnni = 0
        meta.nni = 0
        
        log_sum_weight_subtree = [-Inf]
        valid_subtree = false
        
        pm = 2 * (rand() > 0.5) - 1
        
        rho_bwd = zeros(nl)
        rho_fwd = zeros(nl)
        
        if pm == -1
            transfer!(z_worker, z_minus)
            
            rho_fwd .= rho
            p_fwd_bck .= p_bck_bck
            p_sharp_fwd_bck .= p_sharp_bck_bck

            valid_subtree = buildtree(z_worker,
                z_propose,
                pm,
                depth,
                epsilon,
                logfgrad,
                H0,
                delta,
                nl,
                rho_bwd,
                log_sum_weight_subtree,
                p_sharp_bck_fwd,
                p_sharp_bck_bck,
                p_bck_fwd,
                p_bck_bck,
                meta
            )
        
            transfer!(z_minus, z_worker)

        else

            transfer!(z_worker, z_plus)
            rho_bwd .= rho
            p_bck_fwd .= p_fwd_fwd
            p_sharp_bck_fwd .= p_sharp_fwd_fwd

            valid_subtree = buildtree(z_worker,
                z_propose,
                pm,
                depth,
                epsilon,
                logfgrad,
                H0,
                delta,
                nl,
                rho_fwd,
                log_sum_weight_subtree,
                p_sharp_fwd_bck,
                p_sharp_fwd_fwd,
                p_fwd_bck,
                p_fwd_fwd,
                meta
            )
            
            transfer!(z_plus, z_worker)
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
        log_sum_weight = logsumexp(lsw_t, log_sum_weight)
        nni += meta.nni
        

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

    #@show meta
    v.tune.stepsizeadapter.metro_acc_prob = meta.alpha / meta.nalpha
    #@show nni, tnni
    v.tune.stepsizeadapter.avg_nni = tnni == 0 ? 0.0 : nni/tnni#/(2^depth)#meta.nni / meta.nalpha
    
    
    v.value[1] = z_sample.x
    push!(v.tune.moves, nni)
    push!(v.tune.tree_depth_trace, depth)
    push!(v.tune.acc_p_r, acc_p_r)
    v
end


function buildtree(
    s_worker::Tree_HMC_State,
    s_prob::Tree_HMC_State,
    pm::Int64,
    depth::Integer,
    epsilon::Float64,
    logfgrad::Function,
    H0::Real,
    delta::Float64,
    sz::Int64,
    rho::Vector{Float64},
    log_sum_weight::Vector{Float64},
    p_sharp_beg::Vector{Float64},
    p_sharp_end::Vector{Float64},
    p_beg::Vector{Float64},
    p_end::Vector{Float64},
    meta::PNUTSMeta
    
)

    if depth == 0
        d1 = transfer(s_worker)
        nni = refraction!(d1, pm*epsilon, logfgrad, delta, sz)
        transfer!(s_prob, d1)


        h = hamiltonian(s_worker)
        h = isnan(h) ? -Inf : h
        
        #divergent = H0 < (h + 1000.0)# ? true : false
        divergent = (h-H0) < -1000# ? true : false
        #@show H0, h, divergent#, divergent2
        stat = (h - H0) > 0 ? 1 : exp(h - H0)
        #@show H0, h, stat, divergent, nni
        meta.alpha += stat#H0 - h > 0 ? 1 : exp(H0 - h)
        meta.nalpha += 1
        meta.nni += nni
        if rand() < stat
            meta.accnni += nni
        end
        #@show meta
        log_sum_weight .= logsumexp(log_sum_weight[1], H0-h)
        
        rho .+= s_prob.r
        p_sharp_beg .= s_prob.r
        p_sharp_end .= s_prob.r
        p_beg .=  s_prob.r
        p_end .= s_prob.r
        
        return !divergent
    end

    rho_init = zeros(sz)
    log_sum_weight_init = [-Inf]
    p_init_end = zeros(sz)
    p_sharp_init_end = zeros(sz)
    
    valid_init =
        buildtree(s_worker, s_prob, pm, depth - 1, epsilon, logfgrad, H0, delta, sz, rho_init, log_sum_weight_init,
        p_sharp_beg, p_sharp_init_end, p_beg, p_init_end, meta)
   
    if !valid_init
        return false
    end

    s_final=transfer(s_worker)
    rho_final = zeros(sz)
    log_sum_weight_final = [-Inf]
    p_final_beg = zeros(sz)
    p_sharp_final_beg = zeros(sz)
    
    valid_final =
        buildtree(s_worker, s_final, pm, depth - 1, epsilon, logfgrad, H0, delta, sz, rho_final, log_sum_weight_final,
        p_sharp_final_beg , p_sharp_end, p_final_beg, p_end, meta)
        
    if !valid_final
        return false
    end
       
    log_sum_weight_subtree = logsumexp(log_sum_weight_init[1], log_sum_weight_final[1])
    log_sum_weight[1] = logsumexp(log_sum_weight[1], log_sum_weight_subtree)
    
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


function nouturn(
    rho::Vector{Float64},
    rminus::Vector{Float64},
    rplus::Vector{Float64})::Bool
    return dot(rho, rminus) > 0 && dot(rho, rplus) > 0
end


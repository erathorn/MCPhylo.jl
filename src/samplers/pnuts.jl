#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################
mutable struct PNUTSTune{F<:Function, F2<:Function} <: SamplerTune
    logf::F
    logfgrad::F2
    stepsizeadapter::NUTSstepadapter
    adapt::Bool
    epsilon::Float64
    delta::Float64
    moves::Vector{Int}
    att_moves::Vector{Int}
    tree_depth::Int
    tree_depth_trace::Vector{Int}
    acc_p_r::Vector{Int}
    jitter::Float64



    #PNUTSTune() = new()

    function PNUTSTune(
        x::Vector{T},
        epsilon::Float64,
        logf::F,
        logfgrad::F2;
        target::Real = 0.6,
        tree_depth::Int = 10,
        targetNNI::Float64 = 0.5,
        delta::Float64 = 0.003,
        jitter::Real = 0.0
    ) where {T<:GeneralNode, F, F2}
        @assert 0.0 <= jitter <= 1.0
        new{F, F2}(
            logf,
            logfgrad,
            NUTSstepadapter(
                0,
                0,
                0,
                NUTS_StepParams(0.5, target, 0.05, 0.75, 10, targetNNI)),
            false,
            epsilon,
            delta,
            Int[],
            Int[],
            tree_depth,
            Int[],
            Int[],
            jitter
        )
    end
end

PNUTSTune(
    x::Vector{T},
    logf::Function,
    logfgrad::Function,
    ::NullFunction,
    delta::Float64 = 0.003,
    target::Real = 0.6;
    args...,
) where {T<:GeneralNode} =
    PNUTSTune(x, nutsepsilon(x[1], logfgrad, logf, delta, target), logf, logfgrad; args...)

PNUTSTune(
    x::Vector{T},
    logfgrad::Function,
    delta::Float64,
    target::Real;
    args...,
) where {T<:GeneralNode} =
    PNUTSTune(x, nutsepsilon(x[1], logfgrad, logf, delta, target), logf, logfgrad; args...)

PNUTSTune(x::Vector; epsilon::Real, args...) =
    PNUTSTune(x, epsilon, identity, identity, args...)

const PNUTSVariate = Sampler{PNUTSTune,Vector{T}} where {T<:GeneralNode}


#################### Sampler Constructor ####################


"""
    PNUTS(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for PNUTS sampling. The Parameter is assumed to be
a tree.

Returns a `Sampler{PNUTSTune}` type object.

* params: stochastic node to be updated with the sampler.

* args...: additional keyword arguments to be passed to the PNUTSVariate constructor.
"""
function PNUTS(
    params::ElementOrVector{Symbol};
    delta::Float64 = 0.003,
    target::Float64 = 0.6,
    epsilon::Float64 = -Inf,
    args...,
)

    tune = PNUTSTune(
        GeneralNode[],
        epsilon,
        logpdf!,
        logpdfgrad!;
        delta = delta,
        target = target,
        args...,
    )
    Sampler(params, tune, Symbol[], false)
end


#################### Sampling Functions ####################



#################### Sampling Functions ####################

function sample!(
    v::Sampler{PNUTSTune{F, F2}, Vector{T}},
    logfun::Function;
    grlpdf::Function,
    adapt::Bool = false,
    args...,
)::Sampler{PNUTSTune{F, F2}, Vector{T}} where {T<:GeneralNode, F<:Function, F2<:Function}
    tune = v.tune
    adapter = tune.stepsizeadapter
    
    if adapter.m == 0 && isinf(tune.epsilon)
        tune.epsilon = nutsepsilon(v.value[1], grlpdf, logfun, tune.delta, tune.stepsizeadapter.params.δ)
    end
    setadapt!(v, adapt)
    
    if tune.adapt
        adapter.m += 1
        

        nuts_sub!(v, tune.epsilon, grlpdf, logfun)

        #adapter.metro_acc_prob = adapter.metro_acc_prob# > 1 ? 1 : adapter.metro_acc_prob
        x = dual_averaging(adapter, tune.delta) 
        
        tune.epsilon = exp(x)

    else
        if (adapter.m > 0)
            tune.epsilon = exp(adapter.x_bar)
        end
        nuts_sub!(v, jitter(tune.epsilon, tune.jitter), grlpdf, logfun)
    end
    
    v
end

function dual_averaging(adapter::NUTSstepadapter, delta::F)::F where F<:Real
    const_params = adapter.params
    
    HT = (const_params.δ - adapter.metro_acc_prob) - (const_params.δ_NNI - adapter.NNI_stat)
    
    η = 1.0 / (adapter.m + const_params.t0)
        

    adapter.s_bar = (1.0 - η) * adapter.s_bar + η * HT
    x = const_params.μ - adapter.s_bar * sqrt(adapter.m) / const_params.γ
    # if x < log(0.5 * delta)
    #     adapter.s_bar = ((const_params.μ - log(0.5*delta)) * const_params.γ)/sqrt(adapter.m)
    #     x = log(0.5*delta)
    # end
    

    x_η = adapter.m^-const_params.κ
    adapter.x_bar = (1.0 - x_η) * adapter.x_bar + x_η * x
    return x
end


function setadapt!(v::Sampler{PNUTSTune{F, F2}, Vector{T}}, adapt::Bool)::Sampler{PNUTSTune{F, F2}, Vector{T}} where {F, F2, T}
    tune = v.tune
    if adapt && !tune.adapt
    
        tune.stepsizeadapter.m = 0
        tune.stepsizeadapter.params =
                update_step(tune.stepsizeadapter.params, log(10.0 * tune.epsilon))
    
    end
    tune.adapt = adapt
    v
end


function nuts_sub!(
    v::Sampler{PNUTSTune{F, F2}, Vector{T}},
    epsilon::Float64,
    logfgrad::Function,
    logfun::Function,
)::Sampler{PNUTSTune{F, F2}, Vector{T}} where {F, F2, T}

    #@show epsilon
    x = deepcopy(v.value[1])
    delta = v.tune.delta
    blv = get_branchlength_vector(x)
    r = randn(length(blv))

    blv = get_branchlength_vector(x)
    set_branchlength_vector!(x, molifier.(blv, delta))
    logf, grad = logfgrad(x)
    grad .*= scale_fac.(blv, delta)

    xminus = Tree_HMC_State(deepcopy(x), r[:], grad[:], logf, -1)
    xplus = Tree_HMC_State(deepcopy(x), r[:], grad[:], logf, 1)
    x_prime = Tree_HMC_State(deepcopy(x), r[:], grad[:], logf, 1)
    lu = log(rand())
    logp0 = hamiltonian(xminus)

    nni = 0
    tnni = 0
    j = 0
    n = 1
    sprime = true
    meta = NUTSMeta()
    log_sum_weight = 0.0
    acc_p_r = 0
    ds = 0
    
    while j < v.tune.tree_depth
        pm = rand() > 0.5
        meta.nalpha = 0
        meta.alpha = 0
        meta.l_NNI = 0

        log_sum_weight_subtree = -Inf
        if pm
            
            xprime,_, xminus, nprime, sprime, log_sum_weight_subtree = buildtree(
                xminus,
                xminus.pm,
                j,
                epsilon,
                logfgrad,
                logfun,
                logp0,
                lu,
                delta,
                meta,
                log_sum_weight_subtree,
            )
            
        else
            xprime, xplus,_, nprime, sprime, log_sum_weight_subtree = buildtree(
                xplus,
                xplus.pm,
                j,
                epsilon,
                logfgrad,
                logfun,
                logp0,
                lu,
                delta,
                meta,
                log_sum_weight_subtree,
            )
            
        end
        ds += meta.nalpha
        
        v.tune.stepsizeadapter.metro_acc_prob = meta.alpha / meta.nalpha > 1 ? 1.0 : meta.alpha / meta.nalpha
        #v.tune.stepsizeadapter.NNI_stat = meta.l_NNI/meta.nalpha
        tnni += meta.nni
        if !sprime
            break
        end
        #@show nprime, n
        # sprime is true so checking is not necessary
        if log_sum_weight_subtree > log_sum_weight
            acc_p_r += 1
            v.value[1] = deepcopy(xprime.x)
            nni += meta.nni
        else
            accprob = exp(log_sum_weight_subtree - log_sum_weight)
            if rand() < accprob
                acc_p_r += 1
                v.value[1] = deepcopy(xprime.x)
                nni += meta.nni
            end
        end
        n += nprime
        log_sum_weight = logaddexp(log_sum_weight, log_sum_weight_subtree)
        
        

        j += 1

        s = nouturn(xminus, xplus, epsilon, logfgrad, logfun, delta)

        if !s
            break
        end
    end
    
    v.tune.stepsizeadapter.NNI_stat = ds == 0 ? 0 : meta.att_nni/ds
    push!(v.tune.moves, nni)
    push!(v.tune.att_moves, tnni)
    push!(v.tune.tree_depth_trace, j)
    push!(v.tune.acc_p_r, acc_p_r)
    v
end




function buildtree(
    x::Tree_HMC_State{T},
    pm::Int64,
    j::Integer,
    epsilon::Float64,
    logfgrad::Function,
    logfun::Function,
    logp0::R,
    lu::R,
    delta::Float64,
    meta::NUTSMeta,
    log_sum_weight_subtree::Float64,
)::Tuple{Tree_HMC_State{T},Tree_HMC_State{T},Tree_HMC_State{T},Int,Bool,Float64} where {T<:GeneralNode, R<:Real}


    if j == 0
        
        nni = 0
        att_nni = 0
        if !x.extended
            nni, att_nni = refraction!(x, pm * epsilon, logfgrad, logfun, delta)
        else
            nni = x.nni
            att_nni = x.att_nni
            x.extended = false
        end
        
        logpprime = hamiltonian(x)

        log_sum_weight_subtree = logaddexp(log_sum_weight_subtree, logpprime - logp0)
        
        meta.att_nni += att_nni
        sprime = logp0  < logpprime + 1000.0
        meta.nni += nni
        alphaprime = min(1.0, exp(logpprime - logp0))
        alphaprime = isnan(alphaprime) ? -Inf : alphaprime
        meta.alpha += alphaprime
        nprime = 0
        if (logp0 + lu) < logpprime
            nprime = 1
            meta.accnni += nni
        end
        meta.nalpha += 1
        meta.l_NNI += att_nni
        xprime = transfer(x)
        xplus = transfer(x)
        xminus = transfer(x)
        return xprime, xplus, xminus,nprime, sprime, log_sum_weight_subtree
    else
        log_sum_weight_init = -Inf
        log_sum_weight_final = -Inf
        nprime2 = 0
        
        xprime, xplus, xminus, nprime, sprime, log_sum_weight_init = buildtree(
            x,
            pm,
            j - 1,
            epsilon,
            logfgrad,
            logfun,
            logp0,
            lu,
            delta,
            meta,
            log_sum_weight_init,
        )
        if sprime
            if pm == -1
                worker_final,_,xminus, nprime2, sprime2, log_sum_weight_final = buildtree(
                    x,
                    pm,
                    j - 1,
                    epsilon,
                    logfgrad,
                    logfun,
                    logp0,
                    lu,
                    delta,
                    meta,
                    log_sum_weight_final,
                )
            else
                worker_final,xplus,_, nprime2, sprime2, log_sum_weight_final = buildtree(
                    x,
                    pm,
                    j - 1,
                    epsilon,
                    logfgrad,
                    logfun,
                    logp0,
                    lu,
                    delta,
                    meta,
                    log_sum_weight_final,
                )
            end
        
        
            ls_final = logaddexp(log_sum_weight_init, log_sum_weight_final)
            
            if log_sum_weight_final > ls_final
                transfer!(xprime, worker_final)
            else
                accprob = exp(log_sum_weight_final - ls_final)
                if rand() < accprob
                    transfer!(xprime, worker_final)
                end
            end
            nprime += nprime2
            log_sum_weight_subtree = logaddexp(log_sum_weight_subtree, ls_final)
            
            # overall tree satisfaction
            sprime = sprime2 && nouturn(xplus, xminus, epsilon, logfgrad, logfun, delta)
        end
    end #if j

    return xprime, xplus, xminus, nprime, sprime, log_sum_weight_subtree 
end


function nouturn(
    xminus::Tree_HMC_State{T},
    xplus::Tree_HMC_State{T},
    epsilon::Float64,
    logfgrad::Function,
    logfun::Function,
    delta::Float64,
)::Bool where {T<:GeneralNode}
    _, curr_h = BHV_bounds(xminus.x, xplus.x)
    #curr_h = proj_euc(xminus, xplus)
    if !xminus.extended && !xplus.extended
        temp = Threads.@spawn refraction!(xminus, xminus.pm * epsilon, logfgrad, logfun, delta)
        nni_p, att_nni_p = refraction!(xplus, xplus.pm*epsilon, logfgrad, logfun, delta)
        nni_m, att_nni_m = fetch(temp)
        xminus.extended = true
        xplus.extended = true
        xminus.nni = nni_m
        xplus.nni = nni_p
        xminus.att_nni = att_nni_m
        xplus.att_nni = att_nni_p
    elseif !xminus.extended && xplus.extended
        nni_m, att_nni_m = refraction!(xminus, xminus.pm*epsilon, logfgrad, logfun, delta)
        xminus.extended = true
        xminus.nni = nni_m
        xminus.att_nni = att_nni_m
    elseif xminus.extended && !xplus.extended
        nni_p, att_nni_p = refraction!(xplus, xplus.pm*epsilon, logfgrad, logfun, delta)
        xplus.extended = true
        xplus.nni = nni_p
        xplus.att_nni = att_nni_p
    end
    curr_t_l, _ = BHV_bounds(xminus.x, xplus.x)
    #curr_t_l = proj_euc(xminus, xplus)
    return curr_h <= curr_t_l
end



function proj_euc(xminus, xplus)


    bv1 = MCPhyloTree.get_bipartitions_as_bitvectors(xplus.x)
    bv2 = MCPhyloTree.get_bipartitions_as_bitvectors(xminus.x)
    linds = [n.num for n in get_leaves(xplus.x)]
    blv1 =  get_branchlength_vector(xplus.x)
    blv2 = get_branchlength_vector(xminus.x)
    bl = length(blv1)
    uneq = count(bv1 .!= bv2)
    xdiff=zeros(bl+uneq)
    rminus=zeros(bl+uneq)
    rplus=zeros(bl+uneq)
    ctind = 1
    for i in eachindex(bv1)
        if bv1[i] == bv2[i]
            #if 
            xdiff[i] = i in linds ? abs(blv1[i] - blv2[i]) : blv1[i] - blv2[i]
            rminus[i] = xminus.r[i]
            rplus[i] = xplus.r[i]
        else
            xdiff[i] = blv1[i]
            rplus[i] = xplus.r[i]
            rminus[bl+ctind] = xminus.r[i]
            xdiff[bl+ctind] = -blv2[i]
            ctind +=1
        end
    end
    turbo_dot(xdiff, rplus) >= 0 && turbo_dot(xdiff, rminus) >= 0
end


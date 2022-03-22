#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################
mutable struct PNUTSTune <: SamplerTune
    logf::Union{Function,Missing}
    stepsizeadapter::NUTSstepadapter
    adapt::Bool
    epsilon::Float64
    delta::Float64
    moves::Vector{Int}
    tree_depth::Int
    tree_depth_trace::Vector{Int}
    acc_p_r::Vector{Int}



    PNUTSTune() = new()

    function PNUTSTune(
        x::Vector{T},
        epsilon::Float64,
        logfgrad::Union{Function,Missing};
        target::Real = 0.6,
        tree_depth::Int = 10,
        targetNNI::Float64 = 0.5,
        delta::Float64 = 0.003,
    ) where {T<:GeneralNode}

        new(
            logfgrad,
            NUTSstepadapter(
                0,
                0,
                0,
                NUTS_StepParams(0.5, target, 0.05, 0.75, 10, targetNNI),
            ),
            false,
            epsilon,
            delta,
            Int[],
            tree_depth,
            Int[],
            Int[],
        )
    end
end

PNUTSTune(
    x::Vector{T},
    logfgrad::Function,
    ::NullFunction,
    delta::Float64 = 0.003,
    target::Real = 0.6;
    args...,
) where {T<:GeneralNode} =
    PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta, target), logfgrad; args...)

PNUTSTune(
    x::Vector{T},
    logfgrad::Function,
    delta::Float64,
    target::Real;
    args...,
) where {T<:GeneralNode} =
    PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta, target), logfgrad; args...)

PNUTSTune(x::Vector; epsilon::Real, args...) = PNUTSTune(x, epsilon, missing, args...)

const PNUTSVariate = Sampler{PNUTSTune,T} where {T<:GeneralNode}


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
        logpdfgrad!;
        delta = delta,
        target = target,
        args...,
    )
    Sampler(params, tune, Symbol[], false)
end


#################### Sampling Functions ####################



#################### Sampling Functions ####################

sample!(v::PNUTSVariate; args...) = sample!(v, v.tune.logfgrad; args...)

function sample!(
    v::PNUTSVariate{<:Vector{<:GeneralNode}},
    logfgrad::Function;
    adapt::Bool = false,
    args...,
)
    tune = v.tune
    adapter = tune.stepsizeadapter
    const_params = tune.stepsizeadapter.params

    if adapter.m == 0 && isinf(tune.epsilon)
        tune.epsilon = nutsepsilon(v.value[1], logfgrad, tune.delta, const_params.δ)
    end
    
    setadapt!(v, adapt)
    if tune.adapt
        adapter.m += 1

        nuts_sub!(v, tune.epsilon, logfgrad)

        adaptstat = adapter.metro_acc_prob > 1 ? 1 : adapter.metro_acc_prob

        HT = const_params.δ - adaptstat - (const_params.τ - adapter.avg_nni)

        HT /= 2

        η = 1.0 / (adapter.m + const_params.t0)

        adapter.s_bar = (1.0 - η) * adapter.s_bar + η * HT
        x = const_params.μ - adapter.s_bar * sqrt(adapter.m) / const_params.γ

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


function setadapt!(v::PNUTSVariate{<:Vector{<:GeneralNode}}, adapt::Bool)
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
    v::PNUTSVariate{<:AbstractArray{<:GeneralNode}},
    epsilon::Float64,
    logfgrad::Function,
)
    x = deepcopy(v.value[1])
    delta = v.tune.delta
    blv = get_branchlength_vector(x)
    r = randn(length(blv))

    blv = get_branchlength_vector(x)
    set_branchlength_vector!(x, molifier.(blv, delta))
    logf, grad = logfgrad(x)
    xminus = Tree_HMC_State(deepcopy(x), r, grad, logf)
    xplus = Tree_HMC_State(deepcopy(x), r, grad, logf)

    lu = log(rand())
    logp0 = hamiltonian(xminus)

    nni = 0
    tnni = 0
    j = 0
    n = 1

    meta = NUTSMeta()

    acc_p_r = 0
    while j < v.tune.tree_depth
        pm = 2 * (rand() > 0.5) - 1
        meta.alpha = 0
        meta.nalpha = 0
        meta.accnni = 0
        meta.nni = 0
        if pm == -1

            xminus, _, xprime, nprime, sprime =
                buildtree(xminus, pm, j, epsilon, logfgrad, logp0, lu, delta, meta)

        else

            _, xplus, xprime, nprime, sprime =
                buildtree(xplus, pm, j, epsilon, logfgrad, logp0, lu, delta, meta)


        end#if pm
        v.tune.stepsizeadapter.metro_acc_prob = meta.alpha / meta.nalpha

        tnni += meta.nni
        if !sprime
            break
        end
        # sprime is true so checking is not necessary

        if rand() < nprime / n
            acc_p_r += 1
            v.value[1] = xprime.x
        end

        n += nprime
        nni += meta.nni


        j += 1


        s = nouturn(xminus, xplus, epsilon, logfgrad, delta)

        if !s
            break
        end
    end
    v.tune.stepsizeadapter.avg_nni = tnni == 0 ? 0.0 : nni / tnni
    push!(v.tune.moves, nni)
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
    logp0::Real,
    lu::Real,
    delta::Float64,
    meta::NUTSMeta,
) where {T<:GeneralNode}


    if j == 0
        xprime = transfer(x)
        nni = 0.0
        if !xprime.extended
            nni = refraction!(xprime, pm * epsilon, logfgrad, delta)
        else
            nni = xprime.nni
            xprime.extended = false
        end

        logpprime = hamiltonian(xprime)

        nprime = Int((logp0 + lu) < logpprime)
        sprime = (logp0 + lu) < logpprime + 1000.0

        meta.nni = nni


        xminus = transfer(xprime)
        xplus = transfer(xprime)
        meta.alpha = min(1.0, exp(logpprime - logp0))
        if rand() < meta.alpha
            meta.accnni += nni
        end
        meta.nalpha = 1

    else
        xminus, xplus, xprime, nprime, sprime =
            buildtree(x, pm, j - 1, epsilon, logfgrad, logp0, lu, delta, meta)
        if sprime
            meta1 = NUTSMeta()
            if pm == -1
                xminus, _, xprime2, nprime2, sprime2 =
                    buildtree(xminus, pm, j - 1, epsilon, logfgrad, logp0, lu, delta, meta1)
            else

                _, xplus, xprime2, nprime2, sprime2 =
                    buildtree(xplus, pm, j - 1, epsilon, logfgrad, logp0, lu, delta, meta1)
            end # if pm
            update!(meta, meta1)
            if rand() < nprime2 / (nprime + nprime2)
                transfer!(xprime, xprime2)
            end
            nprime += nprime2

            xminus_bar = transfer(xminus)
            xplus_bar = transfer(xplus)
            sprime = sprime2 && nouturn(xminus_bar, xplus_bar, epsilon, logfgrad, delta)
            if pm == -1
                transfer!(xminus, xminus_bar)
            else
                transfer!(xplus, xplus_bar)
            end

        end #if sprime
    end #if j

    xminus, xplus, xprime, nprime, sprime
end


function nouturn(
    xminus::Tree_HMC_State,
    xplus::Tree_HMC_State,
    epsilon::Float64,
    logfgrad::Function,
    delta::Float64,
)
    curr_l, curr_h = BHV_bounds(xminus.x, xplus.x)

    temp = @spawn refraction!(xminus, -epsilon, logfgrad, delta)
    nni_p = refraction!(xplus, epsilon, logfgrad, delta)
    nni_m = fetch(temp)
    xminus.extended = true
    xplus.extended = true
    xminus.nni = nni_m
    xplus.nni = nni_p
    curr_t_l, curr_t_h = BHV_bounds(xminus.x, xplus.x)
    return curr_h < curr_t_l
end



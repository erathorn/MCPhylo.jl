#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################

mutable struct PNUTSTune <: SamplerTune
    logfgrad::Union{Function,Missing}
    stepsizeadapter::PNUTSstepadapter
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
        targetNNI::Int = 5,
        delta::Float64 = 0.003,
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

const PNUTSVariate = SamplerVariate{PNUTSTune}


#################### Sampler Constructor ####################


"""
    PNUTS(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for PNUTS sampling. The Parameter is assumed to be
a tree.

Returns a `Sampler{PNUTSTune}` type object.

* params: stochastic node to be updated with the sampler.

* args...: additional keyword arguments to be passed to the PNUTSVariate constructor.
"""
function PNUTS(params::ElementOrVector{Symbol}; args...)
    samplerfx = function (model::Model, block::Integer)
        block = SamplingBlock(model, block, true)

        f = let block = block
            (x, sz, ll, gr) -> mlogpdfgrad!(block, x, sz, ll, gr)
        end
        v = SamplerVariate(block, f, NullFunction(); args...)

        sample!(v::PNUTSVariate, f, adapt = model.iter <= model.burnin)

        relist(block, v)
    end
    Sampler(params, samplerfx, PNUTSTune())
end


#################### Sampling Functions ####################

function mlogpdfgrad!(
    block::SamplingBlock,
    x::FNode,
    sz::Int64,
    ll::Bool = false,
    gr::Bool = false,
)::Tuple{Float64,Vector{Float64}}
    grad = Vector{Float64}(undef, sz)
    lp = zero(Float64)

    if gr
        lp, grad = gradlogpdf!(block, x)::Tuple{Float64,Vector{Float64}}
    elseif ll
        lp = logpdf!(block, x)::Float64
    end
    lp, grad
end
sample!(v::PNUTSVariate; args...) = sample!(v, v.tune.logfgrad; args...)

function sample!(v::PNUTSVariate, logfgrad::Function; adapt::Bool = false)
    tune = v.tune
    adapter = tune.stepsizeadapter
    const_params = tune.stepsizeadapter.params
    
    
    setadapt!(v, adapt)
    if tune.adapt
        adapter.m += 1
        

        nuts_sub!(v, tune.epsilon, logfgrad)
        adaptstat = adapter.metro_acc_prob > 1 ? 1 : adapter.metro_acc_prob
        
        HT = const_params.δ - adaptstat
        HT2 = const_params.τ - adapter.avg_nni
        
                
        HT -= abs_adapter(HT2)
        HT /= 2
        
        η = 1.0/(adapter.m + const_params.t0)
        
        adapter.s_bar = (1.0 - η) * adapter.s_bar + η * HT
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


function setadapt!(v::PNUTSVariate, adapt::Bool)
    tune = v.tune
    if adapt && !tune.adapt
        tune.stepsizeadapter.m = 0
        tune.stepsizeadapter.params = update_step(tune.stepsizeadapter.params, log(10.0 * tune.epsilon))
    end
    tune.adapt = adapt
    v
end



function nuts_sub!(v::PNUTSVariate, epsilon::Float64, logfgrad::Function)
    x = deepcopy(v.value[1])
    nl = size(x)[1] - 1
    delta = v.tune.delta
    r = randn(nl)
    #epsilon = 0.01
    
    blv = get_branchlength_vector(x)
    set_branchlength_vector!(x, molifier.(blv, delta))
    logf, grad = logfgrad(x, nl, true, true)
    xminus = Tree_HMC_State(deepcopy(x), r, grad, logf)
    xplus = Tree_HMC_State(deepcopy(x), r, grad, logf)
    
    lu = log(rand())
    logp0 = hamiltonian(xminus)
    

    nni = 0
    j = 0
    n = 1
    
    
    meta = PNUTSMeta()
    directions_plus = rand(Bool, nl)
    directions_minus = rand(Bool, nl)
    acc_p_r = 0
    while j < v.tune.tree_depth
        pm = 2 * (rand() > 0.5) - 1
        
        if pm == -1

            xminus,
            _,
            xprime,
            nprime,
            sprime,_,_ = buildtree(
                xminus,
                pm,
                j,
                epsilon,
                logfgrad,
                logp0,
                lu,
                delta,
                nl,
                meta,
                directions_minus
            )

        else
            
            _,
            xplus,
            xprime,
            nprime,
            sprime,_,_ = buildtree(
                xplus,
                pm,
                j,
                epsilon,
                logfgrad,
                logp0,
                lu,
                delta,
                nl,
                meta,
                directions_plus
            )
            

        end#if pm
        v.tune.stepsizeadapter.metro_acc_prob = meta.alpha / meta.nalpha
        
        v.tune.stepsizeadapter.avg_nni = meta.nni
        
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
        directions_minus = nothing#rand(Bool, nl)
        directions_plus = nothing#rand(Bool, nl)
        s = nouturn(
            xminus,
            xplus,
            epsilon,
            logfgrad,
            delta,
            nl,
            directions_minus,
            directions_plus
        )
        
        if !s
            break
        end
    end

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
    sz::Int64,
    meta::PNUTSMeta,
    directions::Union{Nothing,Vector{Bool}}
) where {T<:FNode}


    if j == 0
        xprime = transfer(x)
        #@show x.r
        nni = refraction!(xprime, pm*epsilon, logfgrad, delta, sz, directions)
        
        logpprime = hamiltonian(xprime)

        nprime = Int((logp0 +lu) < logpprime)
        sprime = (logp0 + lu) < logpprime + 1000.0
        
        #nprime = lu + (logp0 - logpprime) < 0#Int(logu0 < logpprime)

        #sprime = lu + (logp0 - logpprime) < 1000.0
        #sprime = logpprime - (logp0 + lu) > -1000

        xminus = transfer(xprime)
        xplus = transfer(xprime)
        meta.alpha = min(1.0, exp(logpprime - logp0))
        meta.nni = nni
        meta.nalpha = 1
        directions_minus = directions_plus = nothing
    else
        xminus,
        xplus,
        xprime,
        nprime,
        sprime, 
        directions_minus, 
        directions_plus = buildtree(x, pm, j - 1, epsilon, logfgrad, logp0,lu , delta, sz, meta, directions)
        if sprime
            meta1 = PNUTSMeta()
            if pm == -1
                xminus,
                _,
                xprime2,
                nprime2,
                sprime2,
                directions_minus,
                _ = buildtree(
                    xminus,
                    pm,
                    j - 1,
                    epsilon,
                    logfgrad,
                    logp0,
                    lu,
                    delta,
                    sz,
                    meta1,
                    directions_minus
                )
            else
                
                _,
                xplus,
                xprime2,
                nprime2,
                sprime2,
                _,
                directions_plus = buildtree(
                    xplus,
                    pm,
                    j - 1,
                    epsilon,
                    logfgrad,
                    logp0,
                    lu,
                    delta,
                    sz,
                    meta1,
                    directions_plus
                )
            end # if pm
            update!(meta, meta1)
            if rand() < nprime2 / (nprime + nprime2)
                transfer!(xprime, xprime2)
            end
            nprime += nprime2
            #directions_minus = rand(Bool, sz)
            #directions_plus = rand(Bool, sz)
            sprime =
                sprime2 && nouturn(
                    xminus,
                    xplus,
                    epsilon,
                    logfgrad,
                    delta,
                    sz,
                    directions_minus,
                    directions_plus
                )
            
        end #if sprime
    end #if j

    xminus,
    xplus,
    xprime,
    nprime,
    sprime,
    directions_minus,
    directions_plus
    
end


function nouturn(
    xminus::Tree_HMC_State,
    xplus::Tree_HMC_State,
    epsilon::Float64,
    logfgrad::Function,
    delta::Float64,
    sz::Int64,
    directions_minus::Vector{Bool},
    directions_plus::Vector{Bool}
)
    curr_l, curr_h = BHV_bounds(xminus.x, xplus.x)
    xminus_bar = transfer(xminus)
    xplus_bar = transfer(xplus)
     _ = refraction!(
        xminus_bar,
        -epsilon,
        logfgrad,
        delta,
        sz,
        directions_minus
    )
    _ = refraction!(
        xplus_bar,
        epsilon,
        logfgrad,
        delta,
        sz,
        directions_plus
    )

    curr_t_l, curr_t_h = BHV_bounds(xminus_bar.x, xplus_bar.x)
    return curr_h < curr_t_l
end



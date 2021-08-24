#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################

mutable struct PNUTS_Rie_Tune <: SamplerTune
    logfgrad::Union{Function,Missing}
    adapt::Bool
    alpha::Float64
    epsilon::Float64
    epsilonbar::Float64
    gamma::Float64
    Hbar_acc::Float64
    kappa::Float64
    m::Int
    mu::Float64
    nalpha::Int
    t0::Float64
    delta::Float64
    target::Float64
    moves::Vector{Int}
    tree_depth::Int
    nniprime::Float64
    targetNNI::Int
    tree_depth_trace::Vector{Int}


    PNUTS_Rie_Tune() = new()

    function PNUTS_Rie_Tune(
        x::Vector{T},
        epsilon::Float64,
        logfgrad::Union{Function,Missing};
        target::Real = 0.6,
        tree_depth::Int = 10,
        targetNNI::Int = 5,
        delta::Float64 = 0.003,
        jitter::Float64 = 0.0,
    ) where {T<:GeneralNode}

        new(
            logfgrad,
            false,
            0.0,
            epsilon,
            1.0,
            0.05,
            0.0,
            0.75,
            0,
            NaN,
            0,
            10.0,
            delta,
            target,
            Int[],
            tree_depth,
            0,
            targetNNI,
            Int[],
        )
    end
end

PNUTS_Rie_Tune(
    x::Vector{T},
    logfgrad::Function,
    ::NullFunction,
    delta::Float64 = 0.003,
    target::Real = 0.6;
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
    setadapt!(v, adapt)
    if tune.adapt
        tune.m += 1
        tune.nniprime = 0

        nuts_sub!(v, tune.epsilon, logfgrad)
        adaptstat = tune.alpha / tune.nalpha
        adaptstat = adaptstat > 1 ? 1 : adaptstat
        HT = tune.target - adaptstat
        
        HT2 = tune.targetNNI - tune.nniprime / tune.nalpha
        HT -= abs_adapter(HT2)
        HT /= 2
        p = 1.0 / (tune.m + tune.t0)
        tune.Hbar_acc = (1.0 - p) * tune.Hbar_acc + p * HT
        
        tune.epsilon = exp(tune.mu - sqrt(tune.m) * tune.Hbar_acc / tune.gamma)
                
        
        p = tune.m^-tune.kappa
        tune.epsilonbar = exp(p * log(tune.epsilon) + (1.0 - p) * log(tune.epsilonbar))
    else
        if (tune.m > 0)
            tune.epsilon = tune.epsilonbar
        end

        nuts_sub!(v, tune.epsilon, logfgrad)
    end
    v
end


function setadapt!(v::PNUTS_Rie_Variate, adapt::Bool)
    tune = v.tune
    if adapt && !tune.adapt
        tune.m = 0
        tune.mu = log(10.0 * tune.epsilon)
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
    blv = get_branchlength_vector(x)
    set_branchlength_vector!(x, molifier.(blv, delta))
    logf, grad = logfgrad(x, nl, true, true)
    
    #x, r, logf, grad, nni = refraction(x, r, g, epsilon, logfgrad, delta, nl)
    #@show nni
    lu = log(rand())
    logp0 = logf - 0.5 * dot(r)
    logu0 = logp0 + lu#log(rand())
    rminus = rplus = r
    gradminus = gradplus = grad
    rho= r
    xminus = xplus = x
    nni = 0
    j = 0
    n = 1
    s = true
    v.tune.nniprime = 0
    while s && j < v.tune.tree_depth
        
        pm = 2 * (rand() > 0.5) - 1
        rho_bwd = rho_fwd = zeros(nl)
        
        if pm == -1
            

            xminus,
            rminus,
            gradminus,
            _,
            _,
            _,
            xprime,
            nprime,
            sprime,
            alpha,
            nalpha,
            nni1,
            lpp,
            nniprime,
            rho_bwd = buildtree(
                xminus,
                rminus,
                gradminus,
                pm,
                j,
                epsilon,
                logfgrad,
                logp0,
                logu0,
                delta,
                nl,
                lu,
                rho,
            )

        else

            _,
            _,
            _,
            xplus,
            rplus,
            gradplus,
            xprime,
            nprime,
            sprime,
            alpha,
            nalpha,
            nni1,
            lpp,
            nniprime,
            rho_fwd = buildtree(
                xplus,
                rplus,
                gradplus,
                pm,
                j,
                epsilon,
                logfgrad,
                logp0,
                logu0,
                delta,
                nl,
                lu,
                rho
            )

        end#if pm

        v.tune.alpha, v.tune.nalpha = alpha, nalpha
        v.tune.nniprime = nniprime
        
        if !sprime
            break
        end
        # sprime is true so checking is not necessary
        @show nprime/n, nprime, n, j
        if rand() < nprime / n
            println("acc")
            v.value[1] = xprime
        end
        
        n += nprime
        nni += nni1

        rho = rho_bwd .+ rho_fwd
        j += 1
        s = nouturn(
            rho,
            rminus,
            rplus
        )
        @show s
    end

    push!(v.tune.moves, nni)
    push!(v.tune.tree_depth_trace, j)
    v
end


function buildtree(
    x::T,
    r::Vector{Float64},
    grad::Vector{Float64},
    pm::Int64,
    j::Integer,
    epsilon::Float64,
    logfgrad::Function,
    logp0::Real,
    logu0::Real,
    delta::Float64,
    sz::Int64,
    lu::Float64,
    rho::Vector{Float64}
) where {T<:FNode}


    if j == 0

        xprime, rprime, logfprime, gradprime, nni =
            refraction(x, r, grad, pm*epsilon, logfgrad, delta, sz)

        logpprime = logfprime - 0.5 * dot(rprime)
        
        nprime = Int(logu0 < logpprime)
        sprime = logu0 < logpprime + 1000.0
        #@show lu + (logp0 - logpprime), sprime
        xminus = xplus = xprime
        rminus = rplus = rprime
        gradminus = gradplus = gradprime
        alphaprime = min(1.0, exp(logpprime - logp0))
        nniprime = nni
        nalphaprime = 1
        rhoprime = rho + rprime

    else
        rho_init = zeros(sz)

        xminus,
        rminus,
        gradminus,
        xplus,
        rplus,
        gradplus,
        xprime,
        nprime,
        sprime,
        alphaprime,
        nalphaprime,
        nni,
        logpprime,
        nniprime,
        rhoprime =
            buildtree(x, r, grad, pm, j - 1, epsilon, logfgrad, logp0, logu0, delta, sz, lu, rho_init)
        if sprime

            if pm == -1

                xminus,
                rminus,
                gradminus,
                _,
                _,
                _,
                xprime2,
                nprime2,
                sprime2,
                alphaprime2,
                nalphaprime2,
                nni,
                logpprime,
                nniprime2,
                rhoprime = buildtree(
                    xminus,
                    rminus,
                    gradminus,
                    pm,
                    j - 1,
                    epsilon,
                    logfgrad,
                    logp0,
                    logu0,
                    delta,
                    sz,
                    lu,
                    rhoprime
                )
            else
                _,
                _,
                _,
                xplus,
                rplus,
                gradplus,
                xprime2,
                nprime2,
                sprime2,
                alphaprime2,
                nalphaprime2,
                nni,
                logpprime,
                nniprime2,
                rhoprime = buildtree(
                    xplus,
                    rplus,
                    gradplus,
                    pm,
                    j - 1,
                    epsilon,
                    logfgrad,
                    logp0,
                    logu0,
                    delta,
                    sz,
                    lu,
                    rhoprime
                )
            end # if pm

            if rand() < nprime2 / (nprime + nprime2)
                xprime = xprime2
            end
            nprime += nprime2
            sprime = sprime2 && nouturn(rhoprime, rminus, rplus)
            alphaprime += alphaprime2
            nalphaprime += nalphaprime2
            nniprime += nniprime2
        end #if sprime
    end #if j

    xminus,
    rminus,
    gradminus,
    xplus,
    rplus,
    gradplus,
    xprime,
    nprime,
    sprime,
    alphaprime,
    nalphaprime,
    nni,
    logpprime,
    nniprime,
    rhoprime
end


function nouturn(
    rho::Vector{Float64},
    rminus::Vector{Float64},
    rplus::Vector{Float64})::Bool
    return dot(rho, rminus) > 0 && dot(rho, rplus) > 0
end




# #################### Auxilliary Functions ####################

# function nutsepsilon(x::FNode, logfgrad::Function, delta::Float64, target::Float64)

#     x0 = deepcopy(x)
#     n = size(x)[1] - 1

#     # molifier is necessary!
#     blv = get_branchlength_vector(x)
#     set_branchlength_vector!(x, molifier.(blv, delta))
    
#     logf0, gr = logfgrad(x, n, true, true)
    
#     r0 = randn(n)
#     epsilon = 1.0
#     _, rprime, logfprime, _, _ = refraction(x0, r0, gr, epsilon, logfgrad, delta, n)
#     Hp = logfprime - 0.5 * dot(rprime)
#     H0 = logf0 - 0.5 * dot(r0)
#     prob = Hp - H0
#     direction = prob > target ? 1 : -1

#     while direction == 1 ? prob > target : prob < target
#         epsilon = direction == 1 ? 2 * epsilon : 0.5 * epsilon
#         _, rprime, logfprime, _, _ = refraction(x0, r0, gr, epsilon, logfgrad, delta, n)
#         Hp = logfprime - 0.5 * dot(rprime)
#         prob = Hp - H0
#     end
#     epsilon
# end

# @inline function scale_fac(x::T, delta::T) where {T<:Float64}
#     x < delta ? x/delta : 1.0
# end

# @inline function molifier(x::Float64, delta::Float64)::Float64
#     x >= delta ? x : (x^2 + delta^2) / (2.0 * delta)
# end # function

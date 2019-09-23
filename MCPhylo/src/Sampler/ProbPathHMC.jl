



mutable struct ProbPathHMCTune <: SamplerTune
    n_leap::Float64
    stepsz::Float64
    delta::Float64
    logfgrad::Union{Function, Missing}

    ProbPathHMCTune() = new()

    ProbPathHMCTune(n_leap::Float64, stepsz::Float64, delta::Float64) = new(n_leap::Float64, stepsz::Float64, delta::Float64)

    ProbPathHMCTune(n_leap::Float64, stepsz::Float64, delta::Float64, logfgrad::Union{Missing,Function}) = new(n_leap::Float64, stepsz::Float64, delta::Float64, logfgrad)

    function ProbPathHMCTune(x::Vector, n_leap::Float64, stepsz::Float64, delta::Float64, logf::Union{Function, Missing})
      new(n_leap, stepsz, delta, logf)
    end
end # mutable struct

ProbPathHMCTune(x::Vector,n_leap::Float64, stepsz::Float64, delta::Float64 ; args...) =
  ProbPathHMCTune(x, n_leap, stepsz, delta, missing; args...)

const ProbPathVariate = SamplerVariate{ProbPathHMCTune}


function ProbPathHMCSampler(params, pargs...; dtype::Symbol=:forward)

    samplerfx = function(model::Model, block::Integer)

        block = SamplingBlock(model, block, true)
        v = SamplerVariate(block, pargs...)

        sample!(v, block, x -> logpdf!(block, x), x-> gradf!(block, x))

        relist(block, v)
    end # function samplerfx
    Sampler(params, samplerfx, ProbPathHMCTune(pargs...))
end

function sample!(v::ProbPathVariate, block, logf::Function, gradf::Function)

    x1 = v.value[:][1] # copy the original state, so it state can be restored

    n_leap = v.tune.n_leap
    stepsz = v.tune.stepsz
    delta = v.tune.delta

    params = keys(block.model, :block, block.index)
    targets = keys(block.model, :target, block.index)
    a = setdiff(params, targets)

    @assert length(a) == 1

    #ll = logf(v) # log likelihood of the model, including the prior

    #prior = logpdf(block.model, a, false) # log of the prrior


    #tree = block.model[a[1]]

    #probM = randn(size(post_order(tree))[1])

    #currM = deepcopy(probM)
    #currH = ll.+ 0.5*sum(currM.*currM)
    #blens = get_branchlength_vector(tree)
    #probB = deepcopy(blens)

    #for i in 1:n_leap
#
#        fac = scale_factor(v, delta)
#        molify!(v, delta)
#        probM = probM.-stepsz/2.0 .* ((gradf(v).-logpdf(block.model, a, false)).*fac)
#
#        step_nn_att, step_ref_att = refraction(v, probB, probM, true, logf)
#
#
#        set_branchlength_vector!(v.value[1], probB)
#        probM = probM .- stepsz/2.0 .* ((gradf(v).-logpdf(block.model, a, false)).*fac)
#    end

    #probU::Float64 = logf(v)
    #probH = probU + 0.5 * sum(probM.*probM)

    #ratio = currH - probH

    #if ratio <= min(0, log(rand(Uniform(0,1))))
    #    # not successfull
    #    v.value[1] = x1
    #end
    #set_binary!(v.value[1])

    v
end


function refraction(v::ProbPathVariate, probB::Vector{Float64}, probM::Vector{Float64}, surrogate::Bool, logf::Function)

    stepsz = v.tune.stepsz
    delta = v.tune.delta

    postorder = post_order(v.value[1]) # post order and probB are in the same order

    tmpB = @. probB + stepsz * probM
    ref_time = 0.0
    NNI_attempts = 0.0
    ref_attempts = 0.0

    while minimum(tmpB)<=0
        timelist::Vector{Float64} = tmpB./abs.(probM)
        ref_index::Int64 = argmin(timelist)
        temp=(stepsz-ref_time+timelist[ref_index])
        probB = @. probB + temp * probM
        probM[ref_index] *= -1.0
        ref_attempts += 1.0
        if !(postorder[ref_index].nchild == 0)

            if surrogate

                molify!(v, delta)
                U_before_nni::Float64 = logf(v)

                v_copy = deepcopy(v)
                set_binary!(v.value[1])
                tmp_NNI_made = NNI!(v.value[1], postorder[ref_index])


                if tmp_NNI_made != 0

                    molify!(v_copy, delta)
                    U_after_nni::Float64 = logf(v_copy)
                    delta_U = 2.0*(U_after_nni - U_before_nni)
                    my_v::Float64 = probM[ref_index]^2
                    if my_v >= delta_U
                        probM[ref_index] = sqrt(my_v - delta_U)
                        v = v_copy
                        NNI_attempts += 1
                    end
                end
            else
                NNI_attempts += NNI!(v.value[1], ref_index)
            end
        end
        ref_time = stepsz + timelist[ref_index]
        temp = (stepsz-ref_time)
        tmpB = @. probB + temp * probM
    end
    return NNI_attempts, ref_attempts
end # function


function scale_factor(v::SamplerVariate, delta::Float64)::Float64
    mv = minimum(get_branchlength_vector(v.value[1]))
    if mv > delta
        fac = 1.0
    else
        fac = mv/delta
    end
end # function

function molify!(v::SamplerVariate, delta::Float64)
    blens = [molifier(i, delta) for i in get_branchlength_vector(v.value[1])]
    set_branchlength_vector!(v.value[1], blens)
end


"""
    molifier(x::Float64, delta::Float64)::Float64

documentation
"""
function molifier(x::Float64, delta::Float64)::Float64
    if x >= delta
        return delta
    else
        return 1.0/(2.0/delta) * (x*x+delta*delta)
    end
end # function


function gradf!(block::SamplingBlock, x::AbstractVector{T}) where {T<:Real}
    gradlogpdf!(block, x)

end

function mlogpdf!(block::SamplingBlock, x::AbstractVector{T}) where {T<:Real}
  logpdf!(block, x) # logf is already with prior
end

function gradlogpdf!(m::Model, x::AbstractVector{T}, block::Integer=0,
                      transform::Bool=false; dtype::Symbol=:mytree) where {T<:Real}
  GradiantLog(pre_order(x.value[1]), m.nodes[:mypi])
end

function Stochastic(d::AbstractString, f::Function, monitor::Union{Bool, Vector{Int}}=true)
    value = Node()
    fx, src = modelfxsrc(depfxargs, f)
    s = TreeStochastic(value, :nothing, Int[], fx, src, Symbol[],
                      NullUnivariateDistribution())
    setmonitor!(s, monitor)
end

function setmonitor!(d::TreeStochastic, monitor::Bool)
    d.monitor = [1,2]
    d
end

"""
    setinits!(d::TreeVariate, m::model, x::Array)

documentation
"""
function setinits!(d::TreeStochastic, m::Model, x::Array)
    d.value = x
    d.distr = d.eval(m)
    setmonitor!(d, d.monitor)
end # function

function update!(d::TreeStochastic, m::Model)
    d.distr = d.eval(m)
    d
end

function names(d::TreeStochastic, nodekey::Symbol)
    AbstractString["Tree height", "Tree length"]
end

function unlist(d::TreeStochastic)
    y = tree_height(d.value), tree_length(d.value)
    collect(y)
end

function unlist(s::AbstractStochastic, x::Node, transform::Bool=false)
    s.value
end

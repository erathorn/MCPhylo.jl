



mutable struct ProbPathHMCTune <: SamplerTune
    n_leap::Int64
    epsilon::Float64
    delta::Float64
    differ::Union{Symbol, Missing}

    ProbPathHMCTune() = new()

    ProbPathHMCTune(n_leap::Int64, epsilon::Float64, delta::Float64) = new(n_leap::Int64, epsilon::Float64, delta::Float64)

    ProbPathHMCTune(n_leap::Int64, epsilon::Float64, delta::Float64, differ::Union{Missing,Symbol}) = new(n_leap::Int64, epsilon::Float64, delta::Float64, differ)

    function ProbPathHMCTune(x::Vector, n_leap::Int64, epsilon::Float64, delta::Float64, differ::Union{Symbol, Missing})
      new(n_leap, epsilon, delta, differ)
    end
end # mutable struct

ProbPathHMCTune(x::Vector,n_leap::Int64, epsilon::Float64, delta::Float64 ; args...) =
  ProbPathHMCTune(x, n_leap, epsilon, delta, missing; args...)

const ProbPathVariate = SamplerVariate{ProbPathHMCTune}


function ProbPathHMC(params, pargs...)

    samplerfx = function(model::Model, block::Integer)

        block = SamplingBlock(model, block, false)
        v = SamplerVariate(block, pargs...)

        sample!(v, block, x -> logpdf!(block, x), (x, differ) -> gradf!(block, x, differ))

        relist(block, v)
    end # function samplerfx
    Sampler(params, samplerfx, ProbPathHMCTune(pargs...))
end

function mgradient(var::SamplerVariate)
    gradient(var.distr, var.value)
end


function sample!(v::ProbPathVariate, block, logf::Function, gradf::Function)

    x1 = deepcopy(v.value[:][1]) # copy the original state, so it state can be restored

    n_leap = v.tune.n_leap
    stepsz = v.tune.epsilon
    delta = v.tune.delta

    params = keys(block.model, :block, block.index)
    targets = keys(block.model, :target, block.index)

    a = setdiff(params, targets)
    cd_distr = block.model[a[1]].distr

    @assert length(a) == 1

    blens::Vector{Float64} = get_branchlength_vector(block.model[a[1]]) # without root incident length


    currU::Float64 = logf(v) # log likelihood of the model, including the prior

    probM::Vector{Float64} = rand(Normal(),length(blens))
    currM::Vector{Float64} = deepcopy(probM)
    currH::Float64 = currU + 0.5*sum(currM.*currM)
    probB::Vector{Float64} = deepcopy(blens)

    for i in 1:n_leap
        fac::Float64 = scale_factor(v, delta)
        blens_::Vector{Float64} = molify(probB, delta)
        set_branchlength_vector!(v.value[1], blens_)
        relist(block, v)
        l::Vector{Float64} = gradf(v, v.tune.differ)
        l2::Vector{Float64} = gradient(cd_distr, v.value[1])

        probM =  @. probM - stepsz/2.0 * ((l + l2) * fac)

        v, probB, probM = refraction(v, probB, probM, logf, block)

        fac = scale_factor(v, delta)
        blens_ = molify(probB, delta)
        set_branchlength_vector!(v.value[1], blens_)

        relist(block, v)
        l = gradf(v, v.tune.differ)
        l2 = gradient(cd_distr, v.value[1])
        probM =  @. probM - stepsz/2.0 * ((l + l2) * fac)
    end
    set_branchlength_vector!(v.value[1], probB)
    relist(block, v)
    probU::Float64 = logf(v)
    probH = probU + 0.5 * sum(probM.*probM)

    ratio = currH - probH

    if ratio < min(0, log(rand(Uniform(0,1))))
        # not successfull
        v.value[1] = x1
    end

    set_binary!(v.value[1])

    v
end


function refraction(v::ProbPathVariate, probB::Vector{Float64}, probM::Vector{Float64}, logf::Function, block)

    stepsz = v.tune.epsilon
    delta = v.tune.delta

    postorder = post_order(v.value[1])

    tmpB = @. probB + (stepsz * probM)
    ref_time = 0.0


    while minimum(tmpB)<=0
        timelist::Vector{Float64} = tmpB./abs.(probM)
        ref_index::Int64 = argmin(timelist)
        temp=(stepsz-ref_time+timelist[ref_index])
        probB = @. probB + temp * probM

        blens_ = molify(probB, delta)
        set_branchlength_vector!(v.value[1], blens_)
        relist(block, v)
        probM[ref_index] *= -1.0

        if !(postorder[ref_index].nchild == 0)

            U_before_nni::Float64 = logf(v) # still with molified branch length

            v_copy = deepcopy(v)
            tmp_NNI_made = NNI!(v_copy.value[1], postorder[ref_index])

            if tmp_NNI_made != 0
                relist(block, v_copy)
                U_after_nni::Float64 = logf(v_copy)
                delta_U::Float64 = 2.0*(U_after_nni - U_before_nni)
                my_v::Float64 = probM[ref_index]^2

                if my_v >= delta_U
                    probM[ref_index] = sqrt(my_v - delta_U)
                    v = v_copy
                else
                    relist(block, v)
                end # if my_v
            end #if tmpNNI

        end # postorder
        ref_time = stepsz + timelist[ref_index]

        tmpB = @. probB + (stepsz-ref_time) * probM
    end# while
    return v, tmpB, probM
end # function


function scale_factor(v::SamplerVariate, delta::Float64)::Float64
    mv::Float64 = minimum(get_branchlength_vector(v.value[1]))
    fac::Float64 = 0.0
    if mv > delta
        fac = 1.0
    else
        fac = mv/delta
    end
    fac
end # function

function molify(v::Vector{Float64}, delta::Float64)
    return molifier.(v, delta)
end

"""
    molifier(x::Float64, delta::Float64)::Float64

documentation
"""
@inline function molifier(x::Float64, delta::Float64)::Float64
    x >= delta ? x : (x^2+delta^2)/(2.0*delta)
end # function


function gradf!(block::SamplingBlock, x::S, dtype::Symbol=:forward) where {T<:Real, S<:Node}
    #println("here")
    gradlogpdf!(block, x, dtype)

end

function mlogpdf!(block::SamplingBlock, x::AbstractVector{T}) where {T<:Real}
  m = block.model
  index = block.index
  transform = block.transform
  params = keys(m, :block, index)
  targets = keys(m, :target, index)
  m[params] = relist(m, x, params, transform)
  lp = -logpdf(m, setdiff(params, targets), transform)
  for key in targets
    isfinite(lp) || break
    node = m[key]
    update!(node, m)
    lp -= key in params ? logpdf(node, transform) : logpdf(node)
  end
  lp

end


function gradlogpdf(s::AbstractStochastic, x::AbstractArray)
  gradlogpdf_sub(s.distr, x)
end

"""
    setinits!(d::TreeVariate, m::model, x::Array)

documentation
"""

function update!(d::TreeStochastic, m::Model)
    d.distr = d.eval(m)
    d
end

function names(d::TreeStochastic, nodekey::Symbol)
    n_names = [n.num for n in post_order(d.value) if n.root !== true]
    sort!(n_names)
    n_names = vec(AbstractString["node "*string(n) for n in n_names])
    AbstractString["Tree height", "Tree length"]
    vcat(AbstractString["Tree height", "Tree length"], n_names)
end

function names(d::TreeLogical, nodekey::Symbol)
    n_names = [n.num for n in post_order(d.value) if n.root !== true]
    sort!(n_names)
    n_names = vec(AbstractString["node "*string(n) for n in n_names])
    AbstractString["Tree height", "Tree length"]
    vcat(AbstractString["Tree height", "Tree length"], n_names)
end


function unlist(d::TreeStochastic)
    y = tree_height(d.value)
    x = vec([n.height for n in post_order(d.value) if n.root !== true])
    vcat(y, tree_length(d.value), x)

end

function unlist(d::TreeLogical)
    y = tree_height(d.value)
    x = vec([n.height for n in post_order(d.value) if n.root !== true])
    vcat(y, tree_length(d.value), x)
end


function unlist(s::AbstractStochastic, x::Node, transform::Bool=false)
    s.value
end

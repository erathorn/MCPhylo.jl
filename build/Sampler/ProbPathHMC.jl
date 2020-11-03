



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


function molify(v::Vector{Float64}, delta::Float64)
    return molifier.(v, delta)
end

function gradf!(block::SamplingBlock, x::S, dtype::Symbol=:forward) where {S<:Node}
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






function SamplerVariate{T}(x::Symbol, pargs...;kargs...) where T<:ProbPathHMCTune
    println("here")

end

function ProbPathHMC(tree::Symbol, n_leap::Int64, stepsz::Float64)
    println("I am here")
    ProbPathHMCSampler(tree, n_leap, stepsz)
end

function ProbPathHMCSampler(params, pargs...; dtype::Symbol=:forward)

    samplerfx = function(model::Model, block::Integer)
        block = SamplingBlock(model, block, true)
        #v = SamplerVariateP(block, pargs...)
    end # function samplerfx
    Sampler(params, samplerfx, ProbPathHMCTune())
end


function sample!(v::SamplerVariateP)
    tune::ProbPathHMCTune = v.tune # this is the ProbPathHMCTune info
    tree::Node = v.value
    blens_v = get_branchlength_vector(tree)
    blens_copy = blens_v[:]
    curr_U = FelsensteinFunction(post_order(tree), )
    probM = randn(length(blens_copy))
    currM = probM
    curH = curr_U.+ 0.5*sum(curM.*curM)
    probB = blens_copy

    NNI_attempts = 0
    ref_attempts = 0

    for i in 1:tune.n_leap
        tune.stepsz/2.0
        probM = probM.-GradiantLog(pre_order(v), v.pi_).*tune.stepsz/2.0


    end





end # function


function refraction(v, probB, probM)
    tmpB = probB.+ probM.*v.tune.stepsz
    ref_time::Flaot64 = 0.0
    NNI_attempts = 0
    ref_attempts = 0

    while min(tmpB) <= 0
        timelist = tmpB/abs(probM)
        ref_index = argmin(timelist)
        probB = probB .+ (v.tune.stepsz-ref_time+timelist[ref_index]).*probM
        probM[ref_index] *= -1
        ref_attempts += 1

    end

end






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




function sample!(tree::Array{Float64}, data_array::Array{Float64,3}, n_leap::Int64, stepsz::Float64, mypi::Float64, n_c::Int64)

    curr_U = FelsensteinFunction_vec(tree, data_array, mypi, rates, n_c)
    probM = randn(size(tree)[1])
    currM = probM
    curH = curr_U.+ 0.5*sum(curM.*curM)
    probB = blens_copy

    NNI_attempts = 0
    ref_attempts = 0

    tree_c::Array{Float64,2}=deepcopy(tree)

    for i in 1:tune.n_leap
        tune.stepsz/2.0
        probM = probM.-GradiantLog(pre_order(v), v.pi_).*tune.stepsz/2.0
        refraction(tree_c, )

    end





end # function

"""
    refraction(tree::Array{Float64,2}, prob, probm)

documentation
"""
function refraction(tree::Array{Float64,2}, probB::Vector{Float64}, probM::Vector{Float64}, stepsz::Float64, surrogate::bool)
    leaves::Vector{Int64} = get_leaves(tree)
    tmpB = probB .+ stepsz.*probM
    ref_time = 0
    NNI_attempts = 0
    ref_attempts = 0
    while minimum(tmpB)<=0
        timelist::Vector{Float64} = tmpB./abs.(x)
        ref_index::Int64 = argmin(timelist)
        probB = probB .+ (stepsz-ref_time+timelist[ref_index]).*probM
        probM[ref_index] *= -1
        ref_attempts += 1
        if !(ref_index in leaves)
            # do something
            if surrogate is true
                U_before_nni = LogPost(...)

                tree_copy = deepcopy(tree)
                tmp_NNI_made = NNI!(tree, ref_index)

                if tmp_NNI_made != 0
                    U_after_NNI = LogPost(...)
                    delta_U = U_after_NNI - U_before_nni

                    if probM[ref_index]^2 >=2*delta_U

                        probM[ref_index] = sqrt(probM[ref_index]^2 - 2*delta_U)
                        tree = tree_copy
                        NNI_attempts += 1
                    end
                end
            else
                NNI_attempts += NNI!(tree, ref_index)
            end
        end
        ref_time = stepsz + timelist[ref_index]
        tmpB = probB .+ (stepsz-ref_time) .* probM
    end
    return NNI_attempts, ref_attempts
end # function

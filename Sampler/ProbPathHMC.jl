



mutable struct ProbPathHMCTune <: SamplerTune
    n_leap::Int64
    stepsz::Float64

    ProbPathHMCTune() = new()

    ProbPathHMCTune(n_leap::Int64, stpesz::Float64) = new(n_leap::Int64, stpesz::Float64)
end # mutable struct



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


function sample!(tree::Array{Float64}, data_array::Array{Float64,3}, n_leap::Int64, stepsz::Float64, mypi::Float64, n_c::Int64, priordist::T) where {T<:ContinuousMultivariateDistribution}

    curr_U = LogPost(tree, false, delta, rates, n_c, data_array, mypi, priordist)

    probM = randn(size(tree)[1])
    currM = probM
    curH = curr_U.+ 0.5*sum(curM.*curM)
    probB = blens_copy

    NNI_attempts = 0
    ref_attempts = 0

    tree_c::Array{Float64,2}=deepcopy(tree)

    for i in 1:tune.n_leap
        tune.stepsz/2.0
        probM = probM.-GradiantLog(tree_c, data_array, mypi).-logpdf(priordist, probB)
        step_nn_att, step_ref_att = refraction(tree_c, probB, probM, stepsz, surrogate)
        NNI_attempts += step_nn_att
        ref_attempts += step_ref_att

        set_branchlength_vector(tree_c, probB)

        probM = probM .- stepsz/2.0 .* (GradiantLog(tree_c, data_array, mypi).-logpdf(priordist, probB))

    end

    probU = LogPost(tree_c, false, delta, rates, n_c, data, mypi, priordist)
    probH = propU + 0.5 * sum(propM.*probM)

    ratio = currH - probH

    if ratio >= min(0, log(randn(Uniform(0,1))))
        # successfull
        return
    else
        # not successfull
        return
    end



end # function

"""
    refraction(tree::Array{Float64,2}, prob, probm)

documentation
"""
function refraction(tree::Array{Float64,2}, probB::Vector{Float64}, probM::Vector{Float64}, stepsz::Float64, surrogate::Bool)
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
            if surrogate == true
                U_before_nni = LogPost(tree, surrogate, delta, rates, n_c, data_array, mypi, priordist)

                tree_copy = deepcopy(tree)
                tmp_NNI_made = NNI!(tree, ref_index)

                if tmp_NNI_made != 0
                    U_after_NNI = LogPost(tree, surrogate, delta, rates, n_c, data_array, mypi, priordist)
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


"""
    LogPost(tree::Array{Float64,2}, blens::Vector{Float64}, scale::Float64)

documentation
"""
function LogPost(tree::Array{Float64,2}, surrogate::Bool, delta::Float64,
    rates, n_c, data, mypi, priordist::T) where {T<:ContinuousMultivariateDistribution}
    blens = get_branchlength_vector(tree)
    if surrogate
        blens = [molifier(i, delta) for i in blens]
        tree = set_branchlength_vector(tree, blens)
    end
    FelsensteinFunction(tree, data, mypi, rates, n_c)-logpdf(priordist, blens)
end # function


function my_sample!(tree::Node, n_leap::Float64, stepsz::Float64, mypi::Number, n_c::Int64, priordist::Distribution)
    delta = 0.01
    rates = ones(n_c)
    surrogate=true
    curr_U = LogPost(tree, false, delta, rates, n_c, mypi, priordist)

    probM = randn(size(tree)[1])
    currM = deepcopy(probM)
    currH = curr_U.+ 0.5*sum(currM.*currM)
    blens = get_branchlength_vector(tree)
    probB = deepcopy(blens)

    NNI_attempts = 0.0
    ref_attempts = 0.0

    tree_c::Array{Float64,2}=deepcopy(tree)

    for i in 1:n_leap

        probM = probM.-stepsz/2.0 .* GradiantLog(tree_c, mypi, n_c).-_logpdf(priordist, tree_c)
        step_nn_att, step_ref_att = refraction(tree_c, probB, probM, stepsz, surrogate, n_c, delta, rates, priordist, mypi)
        NNI_attempts += step_nn_att
        ref_attempts += step_ref_att

        set_branchlength_vector(tree_c, probB)

        probM = probM .- stepsz/2.0 .* (GradiantLog(tree_c, mypi, n_c).-_logpdf(priordist, tree_c))

    end

    probU::Float64 = LogPost(tree, false, delta, rates, n_c, mypi, priordist)
    probH = probU + 0.5 * sum(probM.*probM)

    ratio = currH - probH

    if ratio >= min(0, log(rand(Uniform(0,1))))
        # successfull
        return tree_c
    else
        # not successfull
        return tree
    end



end # function

"""
    refraction(tree::Array{Float64,2}, prob, probm)

documentation
"""
function refraction(tree::Node, probB::Vector{Float64}, probM::Vector{Float64}, stepsz::Float64, surrogate::Bool,
    n_c::Int64, delta::Float64, rates::Vector{Float64}, priordist::Distribution, mypi::Number)


    postorder = post_order(tree) # post order and probB are in the same order
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
                U_before_nni::Float64 = LogPost(tree, surrogate, delta, rates, n_c, mypi, priordist)

                tree_copy = deepcopy(tree)
                tmp_NNI_made = NNI!(tree, postorder[ref_index])

                if tmp_NNI_made != 0
                    U_after_NNI::Float64 = LogPost(tree, surrogate, delta, rates, n_c, mypi, priordist)
                    delta_U = 2.0*(U_after_NNI - U_before_nni)
                    my_v::Float64 = probM[ref_index]^2
                    if my_v >= delta_U

                        probM[ref_index] = sqrt(my_v - delta_U)
                        tree = tree_copy
                        NNI_attempts += 1
                    end
                end
            else
                NNI_attempts += NNI!(tree, ref_index)
            end
        end
        ref_time = stepsz + timelist[ref_index]
        temp = (stepsz-ref_time)
        tmpB = @. probB + temp * probM
    end
    return NNI_attempts, ref_attempts
end # function


"""
    LogPost(tree::Array{Float64,2}, blens::Vector{Float64}, scale::Float64)

documentation
"""
function LogPost(tree::Node, surrogate::Bool, delta::Float64, rates::Vector{Float64}, n_c::Int64, mypi::Number, priordist::Distribution)::Float64
    blens = get_branchlength_vector(tree)
    if surrogate
        blens = [molifier(i, delta) for i in blens]
        tree = set_branchlength_vector(tree, blens)
    end
    FelsensteinFunction(tree, mypi, rates, n_c)-_logpdf(priordist, tree)
end # function


function my_sample!(tree::Array, data_array::Array{Float64,3}, n_leap::Int64, stepsz::Float64, mypi::Number, n_c::Int64, priordist)
    delta = 0.01
    rates = ones(n_c)
    surrogate=true
    curr_U = LogPost(tree, false, delta, rates, n_c, data_array, mypi, priordist)

    probM = randn(size(tree)[1])
    currM = deepcopy(probM)
    currH = curr_U.+ 0.5*sum(currM.*currM)
    blens = get_branchlength_vector(tree)
    probB = deepcopy(blens)

    NNI_attempts = 0
    ref_attempts = 0

    tree_c::Array{Float64,2}=deepcopy(tree)

    for i in 1:n_leap
        stepsz/2.0
        probM = probM.-stepsz/2.0 .* GradiantLog(tree_c, data_array, mypi, n_c).-_logpdf(priordist, tree_c)
        step_nn_att, step_ref_att = refraction(tree_c, probB, probM, stepsz, surrogate, n_c, data_array, delta, rates, priordist, mypi)
        NNI_attempts += step_nn_att
        ref_attempts += step_ref_att

        set_branchlength_vector(tree_c, probB)

        probM = probM .- stepsz/2.0 .* (GradiantLog(tree_c, data_array, mypi, n_c).-_logpdf(priordist, tree_c))

    end

    probU = LogPost(tree, false, delta, rates, n_c, data_array, mypi, priordist)
    probH = probU + 0.5 * sum(probM.*probM)

    ratio = currH - probH

    if ratio >= min(0, log(rand(Uniform(0,1))))
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
function refraction(tree::Array{Float64,2}, probB::Vector{Float64}, probM::Vector{Float64}, stepsz::Float64, surrogate::Bool,
    n_c, data_array, delta, rates, priordist, mypi)

    leaves::Vector{Int64} = get_leaves(tree)
    tmpB = probB .+ stepsz.*probM
    ref_time = 0
    NNI_attempts = 0
    ref_attempts = 0
    while minimum(tmpB)<=0
        timelist::Vector{Float64} = tmpB./abs.(probM)
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
function LogPost(tree::DenseArray{T}, surrogate::Bool, delta::Float64,
    rates, n_c, data, mypi, priordist) where {T<:Float64}
    blens = get_branchlength_vector(tree)
    if surrogate
        blens = [molifier(i, delta) for i in blens]
        tree = set_branchlength_vector(tree, blens)
    end
    FelsensteinFunction(tree, data, mypi, rates, n_c)-_logpdf(priordist, tree)
end # function

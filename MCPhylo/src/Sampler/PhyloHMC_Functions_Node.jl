




function my_sample!(tree::Node, n_leap::Float64, stepsz::Float64, mypi::Number, n_c::Number, priordist::Distribution)
    delta = 0.01
    rates = ones(n_c)
    surrogate=true
    curr_U = LogPost(tree, false, delta, rates, n_c, mypi, priordist)

    probM = randn(size(post_order(tree))[1])
    currM = deepcopy(probM)
    currH = curr_U.+ 0.5*sum(currM.*currM)
    blens = get_branchlength_vector(tree)
    probB = deepcopy(blens)

    NNI_attempts = 0.0
    ref_attempts = 0.0

    tree_c::Node=deepcopy(tree)


    for i in 1:n_leap

        probM = probM.-stepsz/2.0 .* GradiantLog(pre_order(tree_c), mypi).-_logpdf(priordist, tree_c)
        step_nn_att, step_ref_att = refraction(tree_c, probB, probM, stepsz, surrogate, n_c, delta, rates, priordist, mypi)
        NNI_attempts += step_nn_att
        ref_attempts += step_ref_att

        set_branchlength_vector!(tree_c, probB)

        probM = probM .- stepsz/2.0 .* (GradiantLog(pre_order(tree_c), mypi).-_logpdf(priordist, tree_c))

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
    LogPost(tree::Array{Float64,2}, blens::Vector{Float64}, scale::Float64)

documentation
"""
function LogPost(tree::Node, surrogate::Bool, delta::Float64, rates::Vector{Float64}, n_c::Int64, mypi::Number, priordist::Distribution)::Float64
    if surrogate
        blens = get_branchlength_vector(tree)
        blens = [molifier(i, delta) for i in blens]
        tree = set_branchlength_vector!(tree, blens)
    end
    FelsensteinFunction(post_order(tree), mypi, rates)-_logpdf(priordist, tree)
end # function

function log_grad(tree::Node, mypi::Number, rates::Vector{Float64}, value::Bool, grad::Bool)
    loglik::Float64 = -Inf
    gradll = Vector{Float64}(undef, length(pre_tree))
    if value
        loglik = FelsensteinFunction(post_order(tree), mypi, rates)
    end
    if grad
        pre_tree = pre_order(tree)
        gradll = GradiantLog(pre_tree, mypi)
    end
    return loglik, gradll
end


function jitter(ϵ::R, jitter::R1)::R where {R<:Real,R1<:Real}
    ϵ * (1.0 + jitter * (2.0 * rand() - 1.0))
end

function nutsepsilon(
    x::GeneralNode,
    logfgrad::Function,
    logf::Function,
    delta::Float64,
    target::Float64,
)
    n = size(x)[1] - 1
    log_target = log(target)

    # molifier is necessary!
    #blv = get_branchlength_vector(x)
    #set_branchlength_vector!(x, molifier.(blv, delta))

    logf0, gr = logfgrad(x)

    r0 = randn(n)
    x0 = Tree_HMC_State(deepcopy(x), r0[:], gr[:], logf0)
    x1 = transfer(x0)
    H0 = -hamiltonian(x0)
    epsilon = 1.0
    delta = 2*epsilon
    _ = refraction!(x0, epsilon, logfgrad, logf, delta)
    Hp = -hamiltonian(x0)

    prob = H0 - Hp
    direction = prob > log_target ? 1 : -1

    while direction == 1 ? prob > log_target : prob < log_target
        epsilon = direction == 1 ? 2 * epsilon : 0.5 * epsilon
        
        x1.r = randn(n)
        H0 = -hamiltonian(x1)
        x2 = transfer(x1)
        delta = 2*epsilon
        _ = refraction!(x2, epsilon, logfgrad, logf, delta)

        Hp = -hamiltonian(x2)

        prob = H0 - Hp
    end

    epsilon
end



function nutsepsilon(x::Vector{<:Real}, logfgrad::Function, target::Float64)


    logf0, gr = logfgrad(x)
    target = log(target)

    r0 = randn(size(x))
    x0 = Array_HMC_State(x, r0, gr, logf0)
    x1 = transfer(x0)
    H0 = hamiltonian(x0)
    epsilon = 1.0
    leapfrog!(x0, epsilon, logfgrad)
    Hp = hamiltonian(x0)

    prob = Hp - H0

    direction = prob > target ? 1 : -1
    restart = false
    while direction == 1 ? prob > target : prob < target
        epsilon = direction == 1 ? 2 * epsilon : 0.5 * epsilon
        if isapprox(epsilon, 0)
            throw("The estimated stepsize is too small. ϵ: $epsilon")
        elseif epsilon > 1e7
            throw("The estimated stepsize is too large. ϵ: $epsilon")
        end
        x2 = transfer(x1)
        leapfrog!(x2, epsilon, logfgrad)

        Hp = hamiltonian(x2)

        prob = Hp - H0

    end

    if restart
        @warn "restarting initialization of ϵ"
        epsilon = nutsepsilon(x, logfgrad, exp(target))
    end
    epsilon
end


function nouturn(
    rho::Vector{Float64},
    rminus::Vector{Float64},
    rplus::Vector{Float64},
)::Bool
    return turbo_dot(rminus, rho) > 0 && turbo_dot(rplus, rho) > 0
end


function leapfrog!(d::Array_HMC_State, epsilon::R, logfgrad::Function) where {R<:Real}

    d.r += (0.5 * epsilon) * d.g
    d.x += epsilon * d.r
    d.lf, d.g = logfgrad(d.x)

    d.r += (0.5 * epsilon) * d.g
end

@inline function scale_fac(x::T, delta::T)::T where {T<:Real}
    x < delta ? x / delta : one(T)
end

@inline function molifier(x::T, delta::Real)::T where {T<:Real}
    x >= delta ? x : (x^2 + delta^2) / (2.0 * delta)
end # function

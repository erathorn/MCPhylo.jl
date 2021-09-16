

function nutsepsilon(x::FNode, logfgrad::Function, delta::Float64, target::Float64)
    
    n = size(x)[1] - 1

    # molifier is necessary!
    blv = get_branchlength_vector(x)
    set_branchlength_vector!(x, molifier.(blv, delta))
    
    logf0, gr = logfgrad(x, n, true, true)
    
    r0 = randn(n)
    x0 = Tree_HMC_State(deepcopy(x), r0, gr, logf0)
    x1 = transfer(x0)
    H0 = hamiltonian(x0)
    epsilon = 1.0
    _ = refraction!(x0, epsilon, logfgrad, delta, n, nothing)
    Hp = hamiltonian(x0)
    
    prob = Hp - H0
    direction = prob > target ? 1 : -1

    while direction == 1 ? prob > target : prob < target
        epsilon = direction == 1 ? 2 * epsilon : 0.5 * epsilon
        
        x2 = transfer(x1)
        _ = refraction!(x2, epsilon, logfgrad, delta, n, nothing)
        
        Hp = hamiltonian(x2)
        
        prob = Hp - H0
    end
    
    epsilon
end



@inline function scale_fac(x::T, delta::T)::T where {T<:Real}
    x < delta ? x/delta : 1.0
end

@inline function molifier(x::T, delta::T)::T where {T<:Real}
    x >= delta ? x : (x^2 + delta^2) / (2.0 * delta)
end # function

@inline function abs_adapter(x::R)::Float64 where R <:Real
    x / (1+abs(x))
end
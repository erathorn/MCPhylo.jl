"""
    NUTS(params::ElementOrVector{Symbol}; variant::Symbol=:classic, epsilon::Real = -Inf, kwargs...)

Construct a `Sampler` object for NUTS sampling, with the algorithm´s step size
parameter adaptively tuned during burn-in iterations. Parameters are assumed
to be continuous, but may be constrained or unconstrained. The `variant` parameter controls
whether the Riemannian variant of the NUTS algorithm (https://arxiv.org/abs/1304.1920) or the
classical variant (https://jmlr.org/papers/v15/hoffman14a.html) is used.

Returns a `Sampler{NUTSTune}` type object.

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `variant` : variant of the NUTS sampler. Options are
    * `:classic` : classic version (default)
    * `:riemann` : Riemannian version.

* `epsilon`: step size which can be presupplied. If the standard (-Inf) is chosen, it will be automatically adapted
"""
function NUTS(params::ElementOrVector{Symbol}; variant::Symbol=:classic, epsilon::Real = -Inf, kwargs...)
    if variant == :classic
        return NUTS_classical(params; epsilon=epsilon, kwargs...)
    elseif variant == :riemann
        return NUTS_Rie(params; epsilon=epsilon, kwargs...)      
    else
        throw(ArgumentError("Unrecognized variant: $variant, allowed variants :classic & :riemann"))
    end
end
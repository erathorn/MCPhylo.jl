#################### Monte Carlo Standard Errors ####################
"""
    mcse(x::Vector{T}, method::Symbol=:imse; args...) where {T<:Real}

Compute Monte Carlo standard errors.

Returns the numeric standard error value.

* `x` : time series of values on which to perform calculations.

* `method` : method used for the calculations. Options are
  * `:bm` : batch means, with optional argument `size::Integer=100` determining the number of sequential values to include in each batch. This method requires that the number of values in `x` is at least 2 times the batch size.

  * `:imse` : initial monotone sequence estimator.

  * `:ipse` : initial positive sequence estimator.

* `args...` : additional arguments for the calculation method.
"""
function mcse(x::Vector{T}, method::Symbol = :imse; args...) where {T<:Real}
    method == :bm ? mcse_bm(x; args...) :
    method == :imse ? mcse_imse(x) :
    method == :ipse ? mcse_ipse(x) : throw(ArgumentError("unsupported mcse method $method"))
end

function mcse_bm(x::Vector{T}; size::Integer = 100) where {T<:Real}
    n = length(x)
    m = div(n, size)
    m >= 2 || throw(
        ArgumentError("iterations are < $(2 * size) and batch size is > $(div(n, 2))"),
    )
    mbar = [mean(x[i*size.+(1:size)]) for i = 0:(m-1)]
    sem(mbar)
end

function mcse_imse(x::Vector{T}) where {T<:Real}
    n = length(x)
    m = div(n - 2, 2)
    ghat = autocov(x, [0, 1])
    Ghat = sum(ghat)
    value = -ghat[1] + 2 * Ghat
    for i = 1:m
        Ghat = min(Ghat, sum(autocov(x, [2 * i, 2 * i + 1])))
        Ghat > 0 || break
        value += 2 * Ghat
    end
    sqrt(value / n)
end

function mcse_ipse(x::Vector{T}) where {T<:Real}
    n = length(x)
    m = div(n - 2, 2)
    ghat = autocov(x, [0, 1])
    value = ghat[1] + 2 * ghat[2]
    for i = 1:m
        Ghat = sum(autocov(x, [2 * i, 2 * i + 1]))
        Ghat > 0 || break
        value += 2 * Ghat
    end
    sqrt(value / n)
end

#################### Geweke Diagnostic ####################
"""
    gewekediag(x::Vector{T}; first::Real=0.1, last::Real=0.5,
                etype=:imse, args...) where {T<:Real}


"""
function gewekediag(x::Vector{T}; first::Real=0.1, last::Real=0.5,
                    etype=:imse, args...) where {T<:Real}
  if !(0.0 < first < 1.0)
    throw(ArgumentError("first is not in (0, 1)"))
  elseif !(0.0 < last < 1.0)
    throw(ArgumentError("last is not in (0, 1)"))
  elseif first + last > 1.0
    throw(ArgumentError("first and last proportions overlap"))
  end
  n = length(x)
  x1 = x[1:round(Int, first * n)]
  x2 = x[round(Int, n - last * n + 1):n]
  z = (mean(x1) - mean(x2)) /
      sqrt(mcse(x1, etype; args...)^2 + mcse(x2, etype; args...)^2)
  [round(z, digits=3), round(1.0 - erf(abs(z) / sqrt(2.0)), digits=4)]
end
"""
    gewekediag(c::AbstractChains; first::Real=0.1, last::Real=0.5,
                etype=:imse, args...)

Compute the convergence diagnostic of Geweke [37] for MCMC sampler output. The diagnostic is designed to asses convergence of posterior means estimated with autocorrelated samples. It computes a normal-based test statistic comparing the sample means in two windows containing proportions of the first and last iterations. Users should ensure that there is sufficient separation between the two windows to assume that their samples are independent. A non-significant test p-value indicates convergence. Significant p-values indicate non-convergence and the possible need to discard initial samples as a burn-in sequence or to simulate additional samples.

Returns a `ChainSummary` type object with parameters contained in the rows of the `value` field, and test Z-scores and p-values in the first and second columns. Results are chain-specific.

* `c` : sampler output on which to perform calculations.

* `first` : proportion of iterations to include in the first window.

* `last` : proportion of iterations to include in the last window.

* `etype` : method for computing Monte Carlo standard errors. See `mcse()` for options.

* `args...` : additional arguments to be passed to the `etype` method.

* `x` : vector on which to perform calculations.
"""
function gewekediag(c::AbstractChains; first::Real=0.1, last::Real=0.5,
                    etype=:imse, args...)
  _, p, m = size(c.value)
  vals = Array{Float64}(undef, p, 2, m)
  for j in 1:p, k in 1:m
    vals[j, :, k] = gewekediag(c.value[:, j, k], first=first, last=last,
                               etype=etype; args...)
  end
  hdr = header(c) * "\nGeweke Diagnostic:\nFirst Window Fraction = $first\n" *
        "Second Window Fraction = $last\n"
  ChainSummary(vals, c.names, ["Z-score", "p-value"], hdr)
end

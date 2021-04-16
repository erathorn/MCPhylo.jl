#################### Heidelberger and Welch Diagnostic ####################
"""
    heideldiag(x::Vector{T}; alpha::Real=0.05, eps::Real=0.1,
                etype=:imse, start::Integer=1, args...) where {T<:Real}
"""
function heideldiag(x::Vector{T}; alpha::Real=0.05, eps::Real=0.1,
                    etype=:imse, start::Integer=1, args...) where {T<:Real}
  n = length(x)
  delta = trunc(Int, 0.10 * n)
  y = x[trunc(Int, n / 2):end]
  S0 = length(y) * mcse(y, etype; args...)^2
  i, pvalue, converged, ybar = 1, 1.0, false, NaN
  while i < n / 2
    y = x[i:end]
    m = length(y)
    ybar = mean(y)
    B = cumsum(y) - ybar * collect(1:m)
    Bsq = (B .* B) ./ (m * S0)
    I = sum(Bsq) / m
    pvalue = 1.0 - pcramer(I)
    converged = pvalue > alpha
    if converged
      break
    end
    i += delta
  end
  halfwidth = sqrt(2.0) * erfinv(1.0 - alpha) * mcse(y, etype; args...)
  passed = halfwidth / abs(ybar) <= eps
  [i + start - 2, converged, round(pvalue, digits=4), ybar, halfwidth, passed]
end
"""
    heideldiag(c::AbstractChains; alpha::Real=0.05, eps::Real=0.1,
                etype=:imse, args...)

Compute the convergence diagnostic of Heidelberger and Welch for MCMC sampler output. The diagnostic is designed to assess convergence of posterior means estimated with autocorrelated samples and to determine whether a target degree of accuracy is achieved. A stationarity test is performed for convergence assessment by iteratively discarding 10% of the initial samples until the test p-value is non-significant and stationarity is concluded or until 50% have been discarded and stationarity is rejected, whichever occurs first. Then, a halfwidth test is performed by calculating the relative halfwidth of a posterior mean estimation interval as ``z_{1 - \\alpha / 2} \\hat{s} / |\\bar{\\theta}|``; where ``z`` is a standard normal quantile, ``\\hat{s}`` is the Monte Carlo standard error, and ``\\bar{\\theta}`` is the estimated posterior mean. If the relative halfwidth is greater than a target ratio, the test is rejected. Rejection of the stationarity or halfwidth test suggests that additional samples are needed.

Returns a `ChainSummary` type object with parameters contained in the rows of the `value` field, and numbers of burn-in sequences to discard, whether the stationarity tests are passed (1 = yes, 0 = no), their p-values (``p > \\alpha`` implies stationarity), posterior means, halfwidths of their ``(1 - \\alpha) 100\\%`` estimation intervals, and whether the halfwidth tests are passed (1 = yes, 0 = no) in the columns. Results are chain-specific.

* `c` : sampler output on which to perform calculations.

* `alpha` : significance level for evaluations of stationarity tests and calculations of relative estimation interval halfwidths.

* `eps` : target ratio for the relative halfwidths.

* `etype` : method for computing Monte Carlo standard errors. See `mcse()` for options.

* `args...` : additional arguments to be passed to the `etype` method.

* `x` : vector on which to perform calculations.

* `start` : ??

"""
function heideldiag(c::AbstractChains; alpha::Real=0.05, eps::Real=0.1,
                    etype=:imse, args...)
  _, p, m = size(c.value)
  vals = Array{Float64}(undef, p, 6, m)
  for j in 1:p, k in 1:m
    vals[j, :, k] = heideldiag(c.value[:, j, k], alpha=alpha, eps=eps,
                               etype=etype, start=c.range.start; args...)
  end
  hdr = header(c) * "\nHeidelberger and Welch Diagnostic:\n" *
        "Target Halfwidth Ratio = $eps\nAlpha = $alpha\n"
  ChainSummary(vals, c.names, ["Burn-in", "Stationarity", "p-value", "Mean",
                               "Halfwidth", "Test"], hdr)
end

#################### Raftery and Lewis Diagnostic ####################
"""
    rafterydiag(x::Vector{T}; q::Real=0.025, r::Real=0.005,
                      s::Real=0.95, eps::Real=0.001,
                      range::AbstractRange=1:1:length(x)) where {T<:Real}
"""
function rafterydiag(x::Vector{T}; q::Real=0.025, r::Real=0.005,
                      s::Real=0.95, eps::Real=0.001,
                      range::AbstractRange=1:1:length(x)) where {T<:Real}
  nx = length(x)
  phi = sqrt(2.0) * erfinv(s)
  nmin = ceil(Int, q * (1.0 - q) * (phi / r)^2)
  if nmin > nx
    @warn("At least $nmin samples are needed for specified q, r, and s")
    kthin = burnin = total = NaN
  else
    dichot = Int[(x .<= quantile(x, q))...]
    kthin = 0
    bic = 1.0
    local test, ntest
    while bic >= 0.0
      kthin += 1
      test = dichot[1:kthin:nx]
      ntest = length(test)
      temp = test[1:(ntest - 2)] + 2 * test[2:(ntest - 1)] + 4 * test[3:ntest]
      trantest = reshape(counts(temp, 0:7), 2, 2, 2)
      g2 = 0.0
      for i1 in 1:2, i2 in 1:2, i3 in 1:2
        tt = trantest[i1, i2, i3]
        if tt > 0
          fitted = sum(trantest[:, i2, i3]) * sum(trantest[i1, i2, :]) /
                   sum(trantest[:, i2, :])
          g2 += 2.0 * tt * log(tt / fitted)
        end
      end
      bic = g2 - 2.0 * log(ntest - 2.0)
    end
    tranfinal = counts(test[1:(ntest - 1)] + 2 * test[2:ntest], 0:3)
    alpha = tranfinal[3] / (tranfinal[1] + tranfinal[3])
    beta = tranfinal[2] / (tranfinal[2] + tranfinal[4])
    kthin *= step(range)
    m = log(eps * (alpha + beta) / max(alpha, beta)) /
        log(abs(1.0 - alpha - beta))
    burnin = kthin * ceil(m) + range.start - 1
    n = ((2.0 - alpha - beta) * alpha * beta * phi^2) /
        (r^2 * (alpha + beta)^3)
    keep = kthin * ceil(n)
    total = burnin + keep
  end
  [kthin, burnin, total, nmin, total / nmin]
end
"""
    rafterydiag(c::AbstractChains; q::Real=0.025, r::Real=0.005,
                     s::Real=0.95, eps::Real=0.001)

Compute the convergence diagnostic of Raftery and Lewis for MCMC sampler output. The diagnostic is designed to determine the number of autocorrelated samples required to estimate a specified quantile ``\\theta_q``, such that ``\\Pr(\\theta \\le \\theta_q) = q``, within a desired degree of accuracy. In particular, if ``\\hat{\\theta}_q`` is the estimand and ``\\Pr(\\theta \\le \\hat{\\theta}_q) = \\hat{P}_q`` the estimated cumulative probability, then accuracy is specified in terms of r and s, where ``\\Pr(q - r < \\hat{P}_q < q + r) = s``. Thinning may be employed in the calculation of the diagnostic to satisfy its underlying assumptions. However, users may not want to apply the same (or any) thinning when estimating posterior summary statistics because doing so results in a loss of information. Accordingly, sample sizes estimated by the diagnostic tend to be conservative (too large).

Returns a `ChainSummary` type object with parameters contained in the rows of the value field, and thinning intervals employed, numbers of samples to discard as burn-in sequences, total numbers ``(N)`` to burn-in and retain, numbers of independent samples that would be needed ``(Nmin)``, and dependence factors ``(N / Nmin)`` in the columns. Results are chain-specific.

* `c` : sampler output on which to perform calculations.

* `q` : posterior quantile of interest.

* `r` : margin of error for estimated cumulative probabilities.

* `s` : probability for the margin of error.

* `eps` : tolerance within which the probabilities of transitioning from initial to retained iterations are within the equilibrium probabilities for the chain. This argument determines the number of samples to discard as a burn-in sequence and is typically left at its default value.

* `x` : vector on which to perform calculations.

* `range` : ??
"""
function rafterydiag(c::AbstractChains; q::Real=0.025, r::Real=0.005,
                     s::Real=0.95, eps::Real=0.001)
  _, p, m = size(c.value)
  vals = Array{Float64}(undef, p, 5, m)
  for j in 1:p, k in 1:m
    vals[j, :, k] = rafterydiag(c.value[:, j, k], q=q, r=r, s=s, eps=eps,
                                range=c.range)
  end
  hdr = header(c) * "\nRaftery and Lewis Diagnostic:\n" *
        "Quantile (q) = $q\nAccuracy (r) = $r\nProbability (s) = $s\n"
  ChainSummary(vals, c.names, ["Thinning", "Burn-in", "Total", "Nmin",
                               "Dependence Factor"], hdr)
end

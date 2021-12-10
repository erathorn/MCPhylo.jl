#################### Posterior Statistics ####################
"""
    autocor(c::AbstractChains; lags::Vector=[1, 5, 10, 50],
             relative::Bool=true)

Compute lag-k autocorrelations for MCMC sampler output.

Returns a `ChainSummary` type object with model parameters indexed by the first dimension of `value`, lag-autocorrelations by the second, and chains by the third.

* `c` : sampler output on which to perform calculations.

* `lags` : lags at which to compute autocorrelations.

* `relative` : whether the lags are relative to the thinning interval of the output (`true`) or relative to the absolute iteration numbers (`false`).
"""
function autocor(c::AbstractChains; lags::Vector = [1, 5, 10, 50], relative::Bool = true)
    if relative
        lags *= step(c)
    elseif any(lags .% step(c) .!= 0)
        throw(ArgumentError("lags do not correspond to thinning interval"))
    end
    labels = map(x -> "Lag " * string(x), lags)
    vals = mapslices(x -> autocor(x, lags)', c.value, dims = [1, 2])
    ChainSummary(vals, c.names, labels, header(c))
end
"""
    cor(c::AbstractChains)

Compute cross-correlations for MCMC sampler output.

Returns a `ChainSummary` type object with the first and second dimensions of the `value` field indexing the model parameters between which correlations. Results are for all chains combined.

* `c` : sampler output on which to perform calculations.
"""
function cor(c::AbstractChains)
    ChainSummary(cor(combine(c)), c.names, c.names, header(c))
end
"""
    changerate(c::AbstractChains)

Estimate the probability, or rate per iteration, ``\\Pr(\\theta^i \\ne \\theta^{i-1})`` of a state space change for iterations ``i = 2, \\ldots, N`` in MCMC sampler output. Estimation is performed for each parameter univariately as well as for the full parameter vector multivariately. For continuous output generated from samplers, like Metropolis-Hastings, whose algorithms conditionally accept candidate draws, the probability can be viewed as the acceptance rate.

Returns a `ChainSummary` type object with parameters in the rows of the `value` field, and the estimated rates in the column. Results are for all chains combined.

* `c` : sampler output on which to perform calculations.
"""
function changerate(c::AbstractChains)
    n, p, m = size(c.value)
    r = zeros(Float64, p, 1, 1)
    r_mv = 0.0
    delta = Array{Bool}(undef, p)
    for k = 1:m
        prev = c.value[1, :, k]
        for i = 2:n
            for j = 1:p
                x = c.value[i, j, k]
                dx = x != prev[j]
                r[j] += dx
                delta[j] = dx
                prev[j] = x
            end
            r_mv += any(delta)
        end
    end
    vals = round.([r; r_mv] / (m * (n - 1)), digits = 3)
    ChainSummary(vals, [c.names; "Multivariate"], ["Change Rate"], header(c))
end

describe(c::AbstractChains; args...) = describe(stdout, c; args...)
"""
    describe(io::IO, c::AbstractChains;
                  q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975], etype=:bm, args...)

Compute summary statistics for MCMC sampler output.

Returns results from calls to `summarystats(c, etype, args...)` and `quantile(c, q)` are printed for all chains combined, and a value of `nothing` is returned.

* `c` : sampler output on which to perform calculations.

* `q` : probabilities at which to calculate quantiles.

* `etype` : method for computing Monte Carlo standard errors. See `mcse()` for options.

* `args...` : additional arguments to be passed to the `etype` method.
"""
function describe(
    io::IO,
    c::AbstractChains;
    q::Vector = [0.025, 0.25, 0.5, 0.75, 0.975],
    etype = :bm,
    args...,
)
    ps_stats = summarystats(c; etype = etype, args...)
    ps_quantiles = quantile(c, q = q)
    println(io, ps_stats.header)
    print(io, "Empirical Posterior Estimates:\n")
    show(io, ps_stats)
    print(io, "Quantiles:\n")
    show(io, ps_quantiles)
end
"""
    hpd(x::Vector{T}; alpha::Real=0.05) where {T<:Real}

Compute highest posterior density (HPD) intervals of Chen and Shao [16] for MCMC sampler output. HPD intervals have the desirable property of being the smallest intervals that contain a given probability. However, their calculation assumes unimodal marginal posterior distributions, and they are not invariant to transformations of parameters like central (quantile-based) posterior intervals.

Returns a `ChainSummary` type object with parameters contained in the rows of the value field, and lower and upper intervals in the first and second columns. Results are for all chains combined.

* `c` : sampler output on which to perform calculations.

* `alpha` : the `100 * (1 - alpha)`% interval to compute.

* `x` : vector on which to perform calculations.
"""
function hpd(x::Vector{T}; alpha::Real = 0.05) where {T<:Real}
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n-m+1):n]
    _, i = findmin(b - a)

    [a[i], b[i]]
end
"""
    hpd(c::AbstractChains; alpha::Real=0.05)
"""
function hpd(c::AbstractChains; alpha::Real = 0.05)
    pct = first(showoff([100.0 * (1.0 - alpha)]))
    labels = ["$(pct)% Lower", "$(pct)% Upper"]
    vals = permutedims(
        mapslices(x -> hpd(vec(x), alpha = alpha), c.value, dims = [1, 3]),
        [2, 1, 3],
    )
    ChainSummary(vals, c.names, labels, header(c))
end
"""
    quantile(c::AbstractChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])

Compute posterior quantiles for MCMC sampler output.

Returns a `ChainSummary` type object with parameters contained in the rows of the `value` field, and quantiles in the columns. Results are for all chains combined.

* `c` : sampler output on which to perform calculations.

* `q` : probabilities at which to compute quantiles.
"""
function quantile(c::AbstractChains; q::Vector = [0.025, 0.25, 0.5, 0.75, 0.975])
    labels = map(x -> string(100 * x) * "%", q)
    vals =
        permutedims(mapslices(x -> quantile(vec(x), q), c.value, dims = [1, 3]), [2, 1, 3])
    ChainSummary(vals, c.names, labels, header(c))
end
"""
    summarystats(c::AbstractChains; etype=:bm, args...)

Compute posterior summary statistics for MCMC sampler output.

Returns a `ChainSummary` type object with parameters in the rows of the `value` field; and the sample mean, standard deviation, standard error, Monte Carlo standard error, and effective sample size in the columns. Results are for all chains combined.

* `c` : sampler output on which to perform calculations.

* `etype` : method for computing Monte Carlo standard errors. See `mcse()` for options.

* `args...` : additional arguments to be passed to the `etype` method.
"""
function summarystats(c::AbstractChains; etype = :bm, args...)
    f = x -> [mean(x), std(x), sem(x), mcse(vec(x), etype; args...)]
    labels = ["Mean", "SD", "Naive SE", "MCSE", "ESS"]
    vals = permutedims(mapslices(x -> f(x), c.value, dims = [1, 3]), [2, 1, 3])
    stats = [vals min.((vals[:, 2] ./ vals[:, 4]) .^ 2, size(c.value, 1))]
    ChainSummary(stats, c.names, labels, header(c))
end

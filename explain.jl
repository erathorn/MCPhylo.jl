include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Serialization
using Plots
using StatsPlots
using Gadfly

"""
Mamba Julia
Tutorial

"""
data = rand(Normal(0,1), 5000)

my_data=Dict(:data=>data)

model = Model(
     data = Stochastic(1, (μ, σ) -> Normal(μ, σ), false),
        μ = Stochastic(()->Normal(),true),
        σ = Stochastic(()->Exponential(1), true)
)

inits = [Dict(:data => data,
             :μ => randn(),
             :σ => rand()),
        Dict(:data => data,
            :μ => randn(),
            :σ => rand())]

samplers = [NUTS(:μ),
            Slice(:σ, 0.1)]

setsamplers!(model, samplers)

sim = mcmc(model, my_data, inits, 10000, burnin=500,thin=5, chains=2)

theme(:solarized)
plots = MCPhylo.plot(sim, legend=true)

write("sim.jls", sim)
mysim = read("sim.jls", ModelChains)

using StatsBase

function barplot(c::AbstractChains; legend::Bool=false,
                 position::Symbol=:stack, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plots.Plot}(undef, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    S = unique(c.value[:, i, :])
    n = length(S)
    x = repeat(S, 1, nchains)
    y = zeros(n, nchains)
    for j in 1:nchains
      m = countmap(c.value[:, i, j])
      for k in 1:n
        if S[k] in keys(m)
          y[k, j] = m[S[k]] / nrows
        end
      end
    end
    ymax = maximum(position == :stack ? mapslices(sum, y, dims=2) : y)
    # new plot creation block, based on StatsPlots with a GR backend
    plots[i] = StatsPlots.groupedbar(vec(x), vec(y), bar_position=position,
                                     linecolor=:match,
                                     group=repeat(c.chains, inner=[n]),
                                     legendtitle="Chain", xlabel="Value",
                                     ylabel="Density", title=c.names[i],
                                     legend=pos, ylims=(0.0, ymax),
                                     grid=:dash, gridalpha=0.5)
  end
  return plots
end

function discretediag_sub(c::AbstractChains, frac::Real, method::Symbol,
                          nsim::Int, start_iter::Int, step_size::Int)

  num_iters, num_vars, num_chains = size(c.value)

  vals = zeros(Float64, 3 * (num_chains + 1), num_vars)
  plot_vals_stat = zeros(length(start_iter:step_size:num_iters), num_vars)
  plot_vals_pval = zeros(length(start_iter:step_size:num_iters), num_vars)

  ## Between-chain diagnostic
  X = zeros(Int64, num_iters, num_chains)
  for j in 1:length(num_vars)
    X = convert(Array{Int64, 2}, c.value[:,j,:])
    result = MCPhylo.diag_all(X, method, nsim, start_iter, step_size)
    plot_vals_stat[:,j] = result[1, :] ./ result[2, :]
    plot_vals_pval[:,j] = result[3, :]
    vals[1:3, j] = result[:, end]
  end

  ## Within-chain diagnostic
  x = zeros(Int64, num_iters)
  Y = zeros(Int64, num_iters, 2)
  for j in 1:num_vars
    for k in 1:num_chains
      x = convert(Array{Int64, 1}, c.value[:,j,k])

      idx1 = 1:round(Int64, frac * num_iters)
      idx2 = round(Int64, num_iters - frac * num_iters + 1):num_iters
      x1 = x[idx1]
      x2 = x[idx2]
      n_min = min(length(x1), length(x2))
      Y = [x1[1:n_min] x2[(end - n_min + 1):end]]

      vals[(3 + 3 * (k - 1) + 1):(3 + 3 * (k - 1) + 3), j] =
        MCPhylo.diag_all(Y, method, nsim, n_min, step_size)[:, end]
    end
  end
  return (collect(1:num_vars), vals, plot_vals_stat, plot_vals_pval)

end

function discretediagplot(c::AbstractChains; frac::Real=0.3,
                          method::Symbol=:weiss, nsim::Int=1000,
                          start_iter::Int=100, step_size::Int=10000)

  num_iters, num_vars, num_chains = size(c.value)

  valid_methods = [:hangartner, :weiss, :DARBOOT,
                   :MCBOOT, :billingsley, :billingsleyBOOT]
  if !(method in valid_methods)
    methods_str = join([":$f" for f in valid_methods], ", ")
    throw(ArgumentError("method must be one of ", methods_str))
  end

  if !(0.0 < frac < 1.0)
    throw(ArgumentError("frac must be in (0,1)"))
  end

  if (start_iter > num_iters ) || (step_size > num_iters)
    throw(ArgumentError("start_iter, step_size must be less than $num_iters"))
  end

  V, vals, plot_vals_stat, plot_vals_pval =
    discretediag_sub(c, frac, method, nsim, start_iter, step_size)

  p1 = Gadfly.plot(y=vcat([plot_vals_stat[:,j] for j in 1:length(V)]...),
            x=repeat(collect(c.range[start_iter:step_size:num_iters])/1000,
                     outer=[length(V)]),
            Geom.line,
            Guide.xlabel("Iteration (thousands)", orientation=:horizontal),
            Guide.ylabel("stat/df",orientation=:vertical),
            Scale.color_discrete(), Guide.colorkey(title="Variable"),
            color=repeat(c.names[V],
                         inner=[length(start_iter:step_size:num_iters)]))

   p1_new = Plots.plot(repeat(collect(c.range[start_iter:step_size:num_iters])/
                          1000, outer=[length(V)]),
                       vcat([plot_vals_stat[:,j] for j in 1:length(V)]...),
                       seriestype=:line,
                       group=repeat(c.names[V],
                          inner=[length(start_iter:step_size:num_iters)]),
                       xlabel="Iteration (thousands)", ylabel="stat/df",
                       legendtitle="Variable")

   p2 = Gadfly.plot(y=vcat([plot_vals_pval[:,j] for j in 1:length(V)]...),
            x=repeat(collect(c.range[start_iter:step_size:num_iters])/1000,
                     outer=[length(V)]),
            Geom.line,
            Guide.xlabel("Iteration (thousands)", orientation=:horizontal),
            Guide.ylabel("pval",orientation=:vertical),
            Scale.color_discrete(), Guide.colorkey(title="Variable"),
            color=repeat(c.names[V],
                         inner=[length(start_iter:step_size:num_iters)]))

   p2_new = Plots.plot(repeat(collect(c.range[start_iter:step_size:num_iters])/
                          1000, outer=[length(V)]),
                       vcat([plot_vals_pval[:,j] for j in 1:length(V)]...),
                       seriestype=:line,
                       group=repeat(c.names[V],
                          inner=[length(start_iter:step_size:num_iters)]),
                      xlabel="Iteration (thousands)", ylabel="pval",
                      legendtitle="Variable")


  return [p1, p2]
end



n, p = 25, 10
X = randn(n, p)
beta0 = randn(p)
gamma0 = rand(0:1, p)
y = X * (beta0 .* gamma0) + randn(n)

## Log-transformed Posterior(gamma) + Constant
logf = function(gamma::DenseVector)
  logpdf(MvNormal(X * (beta0 .* gamma), 1.0), y)
end

## MCMC Simulation with Binary Metropolised Gibbs
t = 10000


sim1 = Chains(t, p, names = map(i -> "gamma[$i]", 1:p))
sim2 = Chains(t, p, names = map(i -> "gamma[$i]", 1:p))
gamma1 = BMGVariate(zeros(p), logf)
gamma2 = BMGVariate(zeros(p), logf, k=Vector{Int}[[i] for i in 1:p])
for i in 1:t
  sample!(gamma1)
  sample!(gamma2)
  sim1[i, :, 1] = gamma1
  sim2[i, :, 1] = gamma2
end

sim = cat(sim1, sim2, dims=3)
myplots = discretediagplot(sim, step_size=500)
MCPhylo.draw(myplots)
using Distributed
@everywhere using Mamba, LinearAlgebra, SparseArrays
data = [
 36 27 71  8.1 3.34 11.4 81.5 3243  8.8 42.6 11.7  21  15  59  59  921.87
 35 23 72 11.1 3.14 11.0 78.8 4281  3.5 50.7 14.4   8  10  39  57  997.88
 44 29 74 10.4 3.21  9.8 81.6 4260  0.8 39.4 12.4   6   6  33  54  962.35
 47 45 79  6.5 3.41 11.1 77.5 3125 27.1 50.2 20.6  18   8  24  56  982.29
 43 35 77  7.6 3.44  9.6 84.6 6441 24.4 43.7 14.3  43  38 206  55 1071.29
 53 45 80  7.7 3.45 10.2 66.8 3325 38.5 43.1 25.5  30  32  72  54 1030.38
 43 30 74 10.9 3.23 12.1 83.9 4679  3.5 49.2 11.3  21  32  62  56  934.70
 45 30 73  9.3 3.29 10.6 86.0 2140  5.3 40.4 10.5   6   4   4  56  899.53
 36 24 70  9.0 3.31 10.5 83.2 6582  8.1 42.5 12.6  18  12  37  61 1001.90
 36 27 72  9.5 3.36 10.7 79.3 4213  6.7 41.0 13.2  12   7  20  59  912.35
 52 42 79  7.7 3.39  9.6 69.2 2302 22.2 41.3 24.2  18   8  27  56 1017.61
 33 26 76  8.6 3.20 10.9 83.4 6122 16.3 44.9 10.7  88  63 278  58 1024.89
 40 34 77  9.2 3.21 10.2 77.0 4101 13.0 45.7 15.1  26  26 146  57  970.47
 35 28 71  8.8 3.29 11.1 86.3 3042 14.7 44.6 11.4  31  21  64  60  985.95
 37 31 75  8.0 3.26 11.9 78.4 4259 13.1 49.6 13.9  23   9  15  58  958.84
 35 46 85  7.1 3.22 11.8 79.9 1441 14.8 51.2 16.1   1   1   1  54  860.10
 36 30 75  7.5 3.35 11.4 81.9 4029 12.4 44.0 12.0   6   4  16  58  936.23
 15 30 73  8.2 3.15 12.2 84.2 4824  4.7 53.1 12.7  17   8  28  38  871.77
 31 27 74  7.2 3.44 10.8 87.0 4834 15.8 43.5 13.6  52  35 124  59  959.22
 30 24 72  6.5 3.53 10.8 79.5 3694 13.1 33.8 12.4  11   4  11  61  941.18
 31 45 85  7.3 3.22 11.4 80.7 1844 11.5 48.1 18.5   1   1   1  53  891.71
 31 24 72  9.0 3.37 10.9 82.8 3226  5.1 45.2 12.3   5   3  10  61  871.34
 42 40 77  6.1 3.45 10.4 71.8 2269 22.7 41.4 19.5   8   3   5  53  971.12
 43 27 72  9.0 3.25 11.5 87.1 2909  7.2 51.6  9.5   7   3  10  56  887.47
 46 55 84  5.6 3.35 11.4 79.7 2647 21.0 46.9 17.9   6   5   1  59  952.53
 39 29 75  8.7 3.23 11.4 78.6 4412 15.6 46.6 13.2  13   7  33  60  968.67
 35 31 81  9.2 3.10 12.0 78.3 3262 12.6 48.6 13.9   7   4   4  55  919.73
 43 32 74 10.1 3.38  9.5 79.2 3214  2.9 43.7 12.0  11   7  32  54  844.05
 11 53 68  9.2 2.99 12.1 90.6 4700  7.8 48.9 12.3 648 319 130  47  861.83
 30 35 71  8.3 3.37  9.9 77.4 4474 13.1 42.6 17.7  38  37 193  57  989.27
 50 42 82  7.3 3.49 10.4 72.5 3497 36.7 43.3 26.4  15  18  34  59 1006.49
 60 67 82 10.0 2.98 11.5 88.6 4657 13.5 47.3 22.4   3   1   1  60  861.44
 30 20 69  8.8 3.26 11.1 85.4 2934  5.8 44.0  9.4  33  23 125  64  929.15
 25 12 73  9.2 3.28 12.1 83.1 2095  2.0 51.9  9.8  20  11  26  58  857.62
 45 40 80  8.3 3.32 10.1 70.3 2682 21.0 46.1 24.1  17  14  78  56  961.01
 46 30 72 10.2 3.16 11.3 83.2 3327  8.8 45.3 12.2   4   3   8  58  923.23
 54 54 81  7.4 3.36  9.7 72.8 3172 31.4 45.5 24.2  20  17   1  62 1113.16
 42 33 77  9.7 3.03 10.7 83.5 7462 11.3 48.7 12.4  41  26 108  58  994.65
 42 32 76  9.1 3.32 10.5 87.5 6092 17.5 45.3 13.2  29  32 161  54 1015.02
 36 29 72  9.5 3.32 10.6 77.6 3437  8.1 45.5 13.8  45  59 263  56  991.29
 37 38 67 11.3 2.99 12.0 81.5 3387  3.6 50.3 13.5  56  21  44  73  893.99
 42 29 72 10.7 3.19 10.1 79.5 3508  2.2 38.8 15.7   6   4  18  56  938.50
 41 33 77 11.2 3.08  9.6 79.9 4843  2.7 38.6 14.1  11  11  89  54  946.19
 44 39 78  8.2 3.32 11.0 79.9 3768 28.6 49.5 17.5  12   9  48  53 1025.50
 32 25 72 10.9 3.21 11.1 82.5 4355  5.0 46.4 10.8   7   4  18  60  874.28
 34 32 79  9.3 3.23  9.7 76.8 5160 17.2 45.1 15.3  31  15  68  57  953.56
 10 55 70  7.3 3.11 12.1 88.9 3033  5.9 51.0 14.0 144  66  20  61  839.71
 18 48 63  9.2 2.92 12.2 87.7 4253 13.7 51.2 12.0 311 171  86  71  911.70
 13 49 68  7.0 3.36 12.2 90.7 2702  3.0 51.9  9.7 105  32   3  71  790.73
 35 40 64  9.6 3.02 12.2 82.5 3626  5.7 54.3 10.1  20   7  20  72  899.26
 45 28 74 10.6 3.21 11.1 82.6 1883  3.4 41.9 12.3   5   4  20  56  904.16
 38 24 72  9.8 3.34 11.4 78.0 4923  3.8 50.5 11.1   8   5  25  61  950.67
 31 26 73  9.3 3.22 10.7 81.3 3249  9.5 43.9 13.6  11   7  25  59  972.46
 40 23 71 11.3 3.28 10.3 73.8 1671  2.5 47.4 13.5   5   2  11  60  912.20
 41 37 78  6.2 3.25 12.3 89.5 5308 25.9 59.7 10.3  65  28 102  52  967.80
 28 32 81  7.0 3.27 12.1 81.0 3665  7.5 51.6 13.2   4   2   1  54  823.76
 45 33 76  7.7 3.39 11.3 82.2 3152 12.1 47.3 10.9  14  11  42  56 1003.50
 45 24 70 11.8 3.25 11.1 79.8 3678  1.0 44.8 14.0   7   3   8  56  895.70
 42 33 76  9.7 3.22  9.0 76.2 9699  4.8 42.2 14.5   8   8  49  54  911.82
 38 28 72  8.9 3.48 10.7 79.8 3451 11.7 37.5 13.0  14  13  39  58  954.44
]

pollution = Dict{Symbol, Any}(
  :y => data[:, end],
  :X => mapslices(x -> x / sqrt(var(x)), data[:, 1:(end - 1)], dims=1),
  :p => size(data, 2) - 1
)


## Model Specification
model = Model(

  y = Stochastic(1, (mu, sigma2) -> MvNormal(mu, sqrt(sigma2)), false),

  mu = Logical(1, (alpha, X, theta) -> alpha .+ X * theta, false),

  alpha = Stochastic(() -> Normal(0, 1000)),

  theta = Logical(1, (beta, gamma) -> beta .* gamma),

  beta = Stochastic(1,
    p -> UnivariateDistribution[Normal(0, 1000) for i in 1:p],
    false
  ),

  gamma = Stochastic(1, () -> Bernoulli(0.5)),

  sigma2 = Stochastic(() -> InverseGamma(0.0001, 0.0001))

)


## Gibbs Sampler for alpha and beta
Gibbs_alphabeta = Sampler([:alpha, :beta],
  (alpha, beta, sigma2, X, gamma, y) ->
    begin
      alphabeta_distr = [alpha.distr; beta.distr]
      alphabeta_mean = map(mean, alphabeta_distr)
      alphabeta_invcov = spdiagm(0 => map(d -> 1 / var(d), alphabeta_distr))
      M = [ones(length(y))  X * spdiagm(0 => gamma)]
      Sigma = inv(Symmetric(M' * M / sigma2 + alphabeta_invcov))
      mu = Sigma * (M' * y / sigma2 + alphabeta_invcov * alphabeta_mean)
      alphabeta_rand = rand(MvNormal(mu, Sigma))
      Dict(:alpha => alphabeta_rand[1], :beta => alphabeta_rand[2:end])
    end
)

Gibbs_sigma2 = Sampler([:sigma2],
  (mu, sigma2, y) ->
    begin
      a = length(y) / 2.0 + shape(sigma2.distr)
      b = sum(abs2, y - mu) / 2.0 + scale(sigma2.distr)
      rand(InverseGamma(a, b))
    end
)


## Initial Values
y = pollution[:y]
X = pollution[:X]
p = size(X, 2)

inits = [
  Dict(:y => y, :alpha => mean(y), :gamma => rand(0:1, p),
       :beta => inv(X' * X + Matrix{Float64}(I, p, p)) * X' * y, :sigma2 => var(y)),
  Dict(:y => y, :alpha => 1, :gamma => rand(0:1, p),
       :beta => randn(p), :sigma2 => 1),
  Dict(:y => y, :alpha => 17, :gamma => rand(0:1, p),
       :beta => [15, -15, -10, 5, -10, -5, -10, 10, 40, -5, 0, 0, 0, 20, 5],
       :sigma2 => 1)
]


## Sampling Scheme (without gamma)
scheme0 = [Gibbs_alphabeta, Gibbs_sigma2]

## Binary Hamiltonian Monte Carlo
scheme1 = [BHMC(:gamma, (2 * p + 0.5) * pi); scheme0]
setsamplers!(model, scheme1)
sim1 = mcmc(model, pollution, inits, 10000, burnin=1000, thin=2, chains=3)
describe(sim1)
cs = discretediag(sim1[:, :gamma, :])
showall(cs)

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


  return [p1, p2, p1_new, p2_new]
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
myplots[1]
myplots[2]
myplots[3]
myplots[4]

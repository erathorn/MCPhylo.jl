include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Serialization
using Plots
using StatsPlots

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

sim = mcmc(model, my_data, inits, 5000, burnin=500,thin=5, chains=2)

theme(:solarized)
plots = MCPhylo.plotMC(sim, legend=true)
plots = MCPhylo.plotMC(sim, :mixeddensity)
MCPhylo.draw(plots)


write("sim.jls", sim)
mysim = read("sim.jls", ModelChains)

using StatsBase
function mycontourplot(c::AbstractChains; bins::Integer=100, na...)
  nrows, nvars, nchains = size(c.value)
  # new list initialization
  plots = Plots.Plot[]
  offset = 1e4 * eps()
  n = nrows * nchains
  for i in 1:(nvars - 1)
    X = c.value[:, i, :]
    qx = range(minimum(X) - offset, stop=maximum(X) + offset, length=bins + 1)
    mx = map(k -> mean([qx[k], qx[k + 1]]), 1:bins)
    idx = Int[findfirst(k -> qx[k] <= x < qx[k + 1], 1:bins) for x in X]
    for j in (i + 1):nvars
      Y = c.value[:, j, :]
      qy = range(minimum(Y) - offset, stop=maximum(Y) + offset, length=bins + 1)
      my = map(k -> mean([qy[k], qy[k + 1]]), 1:bins)
      idy = Int[findfirst(k -> qy[k] <= y < qy[k + 1], 1:bins) for y in Y]
      density = zeros(bins, bins)
      for k in 1:n
        density[idx[k], idy[k]] += 1.0 / n
      end

      # new plot creation block, based on Plots with a GR backend
      p = Plots.plot(mx, my, density, seriestype=:contour,
                     colorbar_title="Density", xlabel=c.names[i],
                     ylabel=c.names[j])
      push!(plots, p)
    end
  end
  return plots
end

myplots = mycontourplot(mysim)
MCPhylo.draw(myplots)
myplots[1]
myplots[2]
myplots[3]

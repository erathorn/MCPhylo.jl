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

sim = mcmc(model, my_data, inits, 10000, burnin=500,thin=5, chains=2)

plots = MCPhylo.plotMC(sim)
MCPhylo.draw(plots)

write("sim.jls", sim)
mysim = read("sim.jls", ModelChains)





using StatsBase
function mybarplot(c::AbstractChains; legend::Bool=false,
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
    datapoints = hcat(vec(x), vec(y))
    plots[i] = StatsPlots.groupedbar(datapoints, bar_position = position,
                                     linecolor=:match,
                                     group=repeat(c.chains, inner=[n]),
                                     legendtitle="Chain", xlabel = "Value",
                                     ylabel = "Density", title=c.names[i],
                                     legend = pos, ylims=(0.0, ymax),
                                     grid=:dash, gridalpha=0.5)
  end
  return plots
end
theme(:solarized_light)

myplots = mybarplot(mysim)
myplots[1]
myplots[2]
myplots[3]

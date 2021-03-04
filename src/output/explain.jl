include("../MCPhylo.jl")
using .MCPhylo

trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

"""
plot1 = Plots.plot(trees[1])
plot2 = Plots.plot(trees[1], treetype=:fan, msc=:blue, mc=:yellow, lc=:white,
           bg=:black, tipfont=(7, :lightgreen))
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

sim = mcmc(model, my_data, inits, 1000, burnin=100, thin=5, chains=2)

# default "inner" layout puts plots in a row
pv = plot(sim, [:mean])
# "inner" layout can be manipulated, but usually size has to be adjusted as well
pv = plot(sim, [:mean], layout=(3, 1), size=(800,1500))
# throws an error, as it should for contour (when only one variable is selected)
pv = plot(sim, [:contour], vars=["likelihood"])
# gives a warning for contourplot but shows the other ptypes
pv = plot(sim, [:contour, :density, :mean], vars=["likelihood"], fuse=true)
# specific plot variables are passed successfully
pv = plot(sim, [:autocor, :contour, :density, :mean, :trace],
           maxlag=10, bins=20, trim=(0.1, 0.9), legend=true)
# demonstrate the customizable "outer" layout
pv = plot(sim, [:autocor, :bar, :contour, :mixeddensity, :mean, :trace],
           fuse=true, fLayout=(2,2), fsize=(2750, 2500), linecolor=:match)
# barplot works
pv = plot(sim, [:bar], linecolor=:match, legend=:true, filename="blub.pdf")
# use savefig to save as file; no draw function needed
savefig("test.pdf")

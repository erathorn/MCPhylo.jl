using Revise
# using Distributed
# addprocs(3)
# @everywhere using Pkg
using Pkg
# @everywhere Pkg.activate(".")
Pkg.activate(".")
# @everywhere using MCPhylo
using MCPhylo

### Topology Testing ###
con = MCPhylo.generate_constraints(exc=[(["A", "B", "C", "D", "E"],["F", "G"]), (["a", "b"], ["c"])])
MCPhylo.generate_constraints!(con; exc=[(["a", "b", "c"], ["e", "f"])])
MCPhylo.generate_constraints("./topology.txt")
MCPhylo.generate_constraints!(con, "./topology.txt")
### end Topology testing ###


mt, df = make_tree_with_data("./Example.nex"); # load your own nexus file

mt2 = deepcopy(mt)
randomize!(mt2)
my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :df2 => df,
  :nnodes => size(df)[3],
  :nbase => size(df)[1],
  :nsites => size(df)[2],
);

# model setup
model =  Model(
    df = Stochastic(3, (mtree, mypi) ->  PhyloDist(mtree, mypi, [1.0], [1.0], Restriction), false, false),
    df2 = Stochastic(3, (mtree2, mypi) ->  PhyloDist(mtree2, mypi, [1.0], [1.0], Restriction), false, false),
    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0, 1.0, 0.100, 1.0), true),
    mtree2 = Stochastic(Node(), () -> CompoundDirichlet(1.0, 1.0, 0.100, 1.0), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mtree2 => mt,
    :mypi=> rand(Dirichlet(2,1)),
    :df => my_data[:df],
    :df2 => my_data[:df2],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :a => rand(),
    ),
    Dict{Symbol, Union{Any, Real}}(
        :mtree => mt2,
        :mtree2 => mt2,
        :mypi=> rand(Dirichlet(2,1)),
        :df => my_data[:df],
        :df2 => my_data[:df2],
        :nnodes => my_data[:nnodes],
        :nbase => my_data[:nbase],
        :nsites => my_data[:nsites],
        :a => rand()
        )
]

scheme = [MCPhylo.PNUTS(:mtree, target=0.7, targetNNI=1),
          MCPhylo.PNUTS(:mtree2, target=0.7, targetNNI=1),
           SliceSimplex(:mypi),
          ]

params = SimulationParameters(asdsf=true, freq=25)

setsamplers!(model, scheme)
sim = mcmc(model, my_data, inits, 1000, burnin=100, thin=5, chains=2,
           trees=true, params=params)
MCPhylo.plot_asdsf(sim; legend=true, legendtitlefonthalign=:best, background=:blue)

sim2 = mcmc(sim, 1000, trees=true)
MCPhylo.plot_asdsf(sim2)


#=

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

sim = mcmc(model, my_data, inits, 1000, burnin=500, thin=5, chains=2, trees=true)

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

=#

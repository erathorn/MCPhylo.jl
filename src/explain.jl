using Revise
# using Distributed
# addprocs(3)
# @everywhere using Pkg
using Pkg
# @everywhere Pkg.activate(".")
Pkg.activate(".")
# @everywhere using MCPhylo
using MCPhylo

d = MCPhylo.ParseNewick("./Drav_mytrees_1.nwk")[1]
bd = MCPhylo.BirthDeath(0.2 , 0.23, 0.48)
results1 = [-0.1675285204931672, -0.16316273717762989, -0.20595075254745812, 
            -0.1639099396348209, -0.16246603338187088, -0.16180803461909982, 
            -0.16426413314938396, -0.1637550640877138, -0.16330662842021593, 
            -0.16890167060108943, -0.16436304663505957, -0.17631050764310363, 
            -0.16696457116708, -0.16353869893406714, -0.171275685028695, 
            -0.1634320580060858, -0.17375248059827242, -0.1617230014928498, 
            -0.16727646769660678, -0.1795383581559488, -0.1648908674716657, 
            -0.162017438177007, -0.16189054643044348, -0.16137655321930724, 
            -0.1637879812970666, -0.17318311126257713, -0.17017859269532376, 
            -0.16509378452516307, -0.16132031623282303, -0.16256644553694308, 
            -0.1710818896858957, -0.16961929164080858, -0.1613617725987166, 
            -0.16132067700710798, -0.16259386584710497, -0.16627043949229953, 
            -0.1723098802761037, -0.16337667239460352, -0.20116964123438297, 
            -0.21068699255936466, -0.18467242348418697, -0.18904351151533622, 
            -0.19673144532818773, -0.2669014547836617, -0.16751565759550943, 
            -0.19523798057702957, -0.23574223005760422, -0.1642412674419445, 
            -0.1781905253086248, -0.18075593698449768, -0.16444973804225937, 
            -0.16656017442827428, -0.17058652371102714, -0.1749078841143624, 
            -0.17710995851203762, -0.18326080594126587, -0.18566632464840566, 
            -0.2487611496771592, -0.27787008482416814, -0.1818507463392253, 
            -0.28093334286889254, -0.19424641265044187, -0.2988290412886612, 
            -0.3181268408484228, -0.32928984442510545, -0.3431081569646078, 
            -0.3501474470036245, -0.353430669343034, -0.35764400299316124, 
            -0.1688681421407736, -0.17507299184480057, -0.18386214953649654, 
            -0.36336807459082804]
#=
### Topology Testing ###
con = generate_constraints(exc=[(["A", "B", "C", "D", "E"],["F", "G"]), (["a", "b"], String[])])
con = generate_constraints(mono=[["A", "B", "C", "D", "E"], ["F"]])
generate_constraints!(con; exc=[(["a", "b", "c"], ["e"])])
generate_constraints!(con; mono=[["a", "b", "c"], ["e"]])
generate_constraints("./topology.txt")
generate_constraints!(con, "./topology.txt")
### end Topology testing ###
=#

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
    mtree = Stochastic(Node(), () -> TreeDistribution(CompoundDirichlet(1.0, 1.0, 0.100, 1.0)), true),
    mtree2 = Stochastic(Node(), () -> TreeDistribution(CompoundDirichlet(1.0, 1.0, 0.100, 1.0)), true)
     )
# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(
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

trees = MCPhylo.ParseNewick("./Drav_mytrees_1.nwk")

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


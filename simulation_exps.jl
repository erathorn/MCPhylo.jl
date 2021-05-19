using Distributed
addprocs(3)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using MCPhylo
@everywhere using LinearAlgebra
@everywhere include("simulator.jl")

ntips = 50
leave_names = ["SM_$(i)" for i in 1:ntips]
rt = MCPhylo.create_tree_from_leaves(leave_names)
MCPhylo.number_nodes!(rt)
MCPhylo.set_binary!(rt)
eq_f = rand(4)
eq_f /= sum(eq_f)
rates = [rand()]
data = sample_data(rt, eq_f, rates, JC, 1000)

io = open("gs_params.log","w")
write(io, "eq_f\n")
join(io, eq_f, "\t")
write(io, "\n")
write(io, "rates\n")
join(io, rates, "\t")
write(io, "\n")
write(io, "tree_length\n")
write(io, "$(tree_length(rt))")
write(io, "\n")
write(io, "tree_height\n")
write(io, "$(tree_height(rt))")
close(io)

io = open("gs_tree.nwk","w")
join(io, newick(rt), "\t")
write(io, "\n")
close(io)

mt = deepcopy(rt)
randomize!(mt)
mt2 = deepcopy(mt)
randomize!(mt2)
my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => data,
  :nnodes => size(data)[3],
  :nbase => size(data)[1],
  :nsites => size(data)[2],
);



# model setup
model =  Model(
    df = Stochastic(3, (mtree, mypi, rates) ->  PhyloDist(mtree, mypi, rates, [1.0], JC), false, false),
    mypi = Stochastic(1, () -> Dirichlet(4,1)),
    rates = Stochastic(1, () -> Exponential()),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(Dirichlet(4,1)),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :rates => [rand()],
    ),
    Dict{Symbol, Union{Any, Real}}(
        :mtree => mt2,
        :mypi=> rand(Dirichlet(4,1)),
        :df => my_data[:df],
        :nnodes => my_data[:nnodes],
        :nbase => my_data[:nbase],
        :nsites => my_data[:nsites],
        :rates => [rand()]
        )
    ]

scheme = [PNUTS(:mtree, target=0.8, targetNNI=4),
           SliceSimplex(:mypi),
           Slice(:rates, 1.0)
          ]

setsamplers!(model, scheme);

sp = SimulationParameters(burnin=1250, thin=5, chains=2, trees=true, asdsf=true, freq=5)

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 5000, params=sp)

# request more runs
psrf = max_psrf(sim)
@show psrf
while psrf > 1.1 && sim.stats[end] > 0.01
    global sim, psrf#, bi, gd, gd_values, indices
    @show psrf, sim.stats[end]
    sim = mcmc(sim, 25000)
    psrf = max_psrf(sim)
    to_file(sim, "simulation_")
end

io = open("asdsf_vals.log","w")
join(io, sim.stats, "\n")
close(io)

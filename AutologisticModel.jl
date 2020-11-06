include("./src/MCPhylo.jl")
using .MCPhylo
using Random
import Distributions: logpdf, DiscreteMatrixDistribution
include("./AutologistcDistr.jl")
include("./linguistic_features.jl")
include("./neighbour_graphs.jl")

Random.seed!(42)


wo_nmat = create_nmat(wo_dmat, 500)
wo_nmat = Int64.(wo_nmat)
wo_lmat = create_linguistic_nmat(wo_data)

wo_values
spatial_sums = get_neighbour_sums(wo_values, wo_ngraph)
linguistic_sums = get_neighbour_sums(wo_values, wo_lgraph)

# model setup
model =  Model(
    df = Stochastic(3, (Hcon, Hw, Vcon, Vw, Usum, Uw) ->
                    AutologisticDistr(Hcon, Hw, Vcon, Vw, Usum, Uw), false, false),

    Hcon = Logical(2, (graphweights) -> weightGraphandCount(graphweights, OG), false),
    graphweights = Stochastic(1, ()-> Gamma(), true),

    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

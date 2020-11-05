


Hcon = Array{Float64, 2}(nlangs, nfeat)

function weightGraphandCount(Graph, weights)
    nothing
end

# model setup
model =  Model(
    df = Stochastic(3, (Hcon, Hw, Vcon, Vw, Usum, Uw) ->
                    ALD(Hcon, Hw, Vcon, Vw, Usum, Uw), false, false),

    Hcon = Logical(2, (graphweights) -> weightGraphandCount(graphweights, OG), false),
    graphweights = Stochastic(1, ()-> Gamma(), true),

    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

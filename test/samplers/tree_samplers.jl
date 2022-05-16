# load data
Random.seed!(123)

mt, df = make_tree_with_data("./likelihood/simudata.nex", binary=true)

mt2 = deepcopy(mt)

randomize!(mt2)
my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :nnodes => size(df)[3],
  :nbase => size(df)[1],
  :nsites => size(df)[2],
)


# model setup
model =  Model(
    df = Stochastic(3, (mtree, mypi) ->  PhyloDist(mtree, mypi, [1.0], [1.0], JC), false),
    mypi = Stochastic(1, () -> Dirichlet(4,1)),
    mtree = Stochastic(Node(), () -> TreeDistribution(CompoundDirichlet(1.0,1.0,0.100,1.0)), true)
     )
# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(Dirichlet(4,1)),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :a => rand(),
    ),
    Dict{Symbol, Union{Any, Real}}(
        :mtree => mt2,
        :mypi=> rand(Dirichlet(4,1)),
        :df => my_data[:df],
        :nnodes => my_data[:nnodes],
        :nbase => my_data[:nbase],
        :nsites => my_data[:nsites],
        :a => rand()
        )
    ]


@testset "PNUTS" begin
    Random.seed!(123)

    scheme = [PNUTS(:mtree, target=0.8, targetNNI=0.6 tree_depth=5), SliceSimplex(:mypi)]

    setsamplers!(model, scheme)
    sim = mcmc(
        model,
        my_data,
        inits,
        1000,
        burnin = 50,
        thin = 1,
        chains = 2,
        verbose=false,
        trees=true
    )

    mcmc(sim, 1500)

    r_m = summarystats(sim).value[5, 1, 1]
    r_sd = summarystats(sim).value[5, 2, 1]
    @test r_m - r_sd <= 1.8 <= r_m + r_sd 
    leafset_1 = sort!([l.name for l in get_leaves(ParseNewick(sim.trees[50]))]) 
    leafset_2 = sort!([l.name for l in get_leaves(ParseNewick(sim.trees[750]))]) 
    @test leafset_1 == leafset_2 

    @testset "random draw from Abstract Stochastic" begin
        result = rand(sim.model[:mypi], 5)
        @test size(result) == (4, 5) 
        result = rand(sim.model[:mypi])
        @test isa(result, Vector{Float64})
        @test length(result) == 4
    end
end

@testset "RWM_Univariate" begin
    Random.seed!(123)

    scheme = [RWM(:mtree, :all), Slice(:mtree, 1.0, Univariate), SliceSimplex(:mypi)]

    setsamplers!(model, scheme)

    sp = SimulationParameters(burnin=50, chains=2, verbose=false, trees=true, asdsf=true)
    sim = mcmc(
        model,
        my_data,
        inits,
        1000,
        params=sp
    )
    
    r_m = summarystats(sim).value[5, 1, 1]
    r_sd = summarystats(sim).value[5, 2, 1]
    @test r_m - r_sd <= 1.8 <= r_m + r_sd

    @testset "ASDSF on chain & on the fly" begin
        asdsf = ASDSF(sim; freq=50)
        @test sim.stats[:, 1] == asdsf[1]
        @test all(x -> x <= 1, sim.stats)
    end
end


@testset "RWM_Multivariate" begin
    Random.seed!(123)

    scheme = [RWM(:mtree, :all), Slice(:mtree, 1.0, Multivariate), SliceSimplex(:mypi)]

    setsamplers!(model, scheme)

    sp = SimulationParameters(burnin=50, chains=2, verbose=false, trees=true, asdsf=false)
    sim = mcmc(
        model,
        my_data,
        inits,
        2000,
        params=sp
    )
    
    r_m = summarystats(sim).value[5, 1, 1]
    r_sd = summarystats(sim).value[5, 2, 1]
    @test r_m - r_sd <= 1.8 <= r_m + r_sd
end


@testset "PPHMC" begin
    Random.seed!(123)

    scheme = [PPHMC(:mtree, 0.001, 10), SliceSimplex(:mypi)]

    setsamplers!(model, scheme)
        sim = mcmc(
            model,
            my_data,
            inits,
            10000,
            burnin = 2000,
            thin = 1,
            chains = 2,
            verbose=false,
            trees=true
        )
    
    r_m = summarystats(sim).value[5, 1, 1]
    r_sd = summarystats(sim).value[5, 2, 1]
    @test r_m - r_sd <= 1.8 <= r_m + r_sd
end

@testset "abc exponential" begin
    gs_p = randn()^2
    
    data = Dict{Symbol,Any}(
        :y => rand(Exponential(gs_p), 1000)
    )
    
    mod = Model(
        y=Stochastic(1, (l) -> Exponential(l), false),
        l=Stochastic(() -> LogNormal(), true)
    )
    
    inits = [
        Dict{Symbol,Any}(
          :y => data[:y],
          :l => rand(LogNormal()),
            )
        for i in 1:3
    ]

    f = x -> [mean(x), std(x)]
    scheme = [ABC(:l, 0.01, f, 0.01, nsim = 3, maxdraw=15)]

    setsamplers!(mod, scheme)

    sim = mcmc(mod, data, inits, 100000, burnin=25000, thin=10, chains=2, verbose=false)
    
    while maximum(gelmandiag(sim).value[1,2,:]) > 1.1
        @show maximum(gelmandiag(sim).value[1,2,:])
        sim = mcmc(sim, 10000)
    end
    
    m = mean(sim.value[:,1,:])
    sd = std(sim.value[:,1,:])

    @test m - sd <= gs_p <= m + sd
end


@testset "abc Normal" begin
    gs_m = randn()
    gs_s = randn()^2
    data = Dict{Symbol,Any}(
        :y => rand(Normal(gs_m, gs_s), 1000)
    )
    
    mod = Model(
        y=Stochastic(1, (m, s) -> Normal(m, s), false),
        m=Stochastic(() -> Normal(), true),
        s=Stochastic(() -> Exponential(), true)
    )
    
    inits = [
        Dict{Symbol,Any}(
          :y => data[:y],
          :m => gs_m,
          :s => rand(),
        )
        for i in 1:3
    ]
    
    f = x -> [mean(x), std(x)]
    scheme = [ABC(:m, 0.01, f, 0.001, nsim = 3, maxdraw=15),
              ABC(:s, 0.01, f, 0.001, nsim = 3, maxdraw=15)]

    setsamplers!(mod, scheme)

    sim = mcmc(mod, data, inits, 100000, burnin=25000, thin=10, chains=2, verbose=false)
    
    while maximum(gelmandiag(sim).value[1:2,2,:]) > 1.1
        @show maximum(gelmandiag(sim).value[1:2,2,:])
        sim = mcmc(sim, 10000)
    end
    
    m = mean(sim.value[:,1,:])
    sd = std(sim.value[:,1,:])

    @test m - sd <= gs_s <= m + sd

    m = mean(sim.value[:,2,:])
    sd = std(sim.value[:,2,:])
    
    @test m - sd <= gs_m <= m + sd
end
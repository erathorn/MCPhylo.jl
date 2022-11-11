@testset "Sampler Construction" begin

    s = DMH(:x, 15)
    @test typeof(s) <: Sampler
    
    s = DMH(:x, 15, scale = 0.1)
    @test typeof(s) <: Sampler

    s = MISS(:x)
    @test typeof(s) <: Sampler

    s = Empirical(:x, 5)
    @test typeof(s) <: Sampler

    s = ABC(:x, 1.0, mean, 0.01)
    @test typeof(s) <: Sampler

end
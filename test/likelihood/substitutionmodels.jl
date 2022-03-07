@testset "SubstitutionModels.jl" begin
    @testset "Restriction" begin
        target = ([-0.9191450300180579 -0.7071067811865475; 0.3939192985791677 -0.7071067811865476], 
                  [-1.0, 0.0], [-0.7615773105863908 0.7615773105863908; -0.4242640687119285 -0.9899494936611665], 
                  2.380952380952381)
        @test Restriction([0.3, 0.7], [0.0]) == target
        @test_throws AssertionError MCPhylo.Restriction(ones(3) / 3, [0.1])
    end   
end
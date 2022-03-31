@testset "varaite.jl" begin

    @testset "ScalarLogical 1" begin
        x = rand()
        y = rand()
        @test MCPhylo.ScalarLogical(x) + y == x+y
        @test MCPhylo.ScalarLogical(x) - y == x-y
        @test MCPhylo.ScalarLogical(x) * y == x*y
        @test MCPhylo.ScalarLogical(x) / y == x/y
        @test MCPhylo.ScalarLogical(x) \ y == x\y
    end
    @testset "ScalarLogical 2" begin
        x = rand()
        y = rand()
        @test MCPhylo.ScalarLogical(x) + MCPhylo.ScalarLogical(y) == x+y
        @test MCPhylo.ScalarLogical(x) - MCPhylo.ScalarLogical(y) == x-y
        @test MCPhylo.ScalarLogical(x) * MCPhylo.ScalarLogical(y) == x*y
        @test MCPhylo.ScalarLogical(x) / MCPhylo.ScalarLogical(y) == x/y
        @test MCPhylo.ScalarLogical(x) \ MCPhylo.ScalarLogical(y) == x\y
    end

    @testset "Stochastic Variate 1" begin
        x = rand()
        y = rand()
        s1 = Stochastic(() -> Normal(), true)
        s1.value = x
        
        @test s1 + y == x+y
        @test s1 - y == x-y
        @test s1 * y == x*y
        @test s1 / y == x/y
        @test s1 \ y == x\y    
    end
    @testset "Stochastic Variate 2" begin
        x = rand()
        y = rand()
        s1 = Stochastic(() -> Normal(), true)
        s1.value = x
        s2 = Stochastic(() -> Normal(), true)
        s2.value = y
        @test s1 + s2 == x+y
        @test s1 - s2 == x-y
        @test s1 * s2 == x*y
        @test s1 / s2 == x/y
        @test s1 \ s2 == x\y
    end

    
    @testset "Array Variate 1" begin
        x = rand(2)
        y = rand(2)
        s1 = Stochastic(1, () -> Normal(), true)
        s1.value = x
        @test s1 + y == x+y
        @test s1 - y == x-y
        @test s1 / y == x/y
        @test s1 \ y == x\y
    end
    @testset "Array Variate 2" begin
        x = rand(2)
        y = rand(2)
        s1 = Stochastic(1, () -> Normal(), true)
        s1.value = x
        s2 = Stochastic(1, () -> Normal(), true)
        s2.value = y
        
        @test s1 + s2 == x+y
        @test s1 - s2 == x-y
        @test s1 / s2 == x/y
        @test s1 \ s2 == x\y
    end
end
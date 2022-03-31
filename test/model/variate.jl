@testset "variate.jl" begin

    @testset "conversion scalar" begin
        y = Logical(()->1, true)
        y.value = 1
        @test convert(Bool, y) == true
        @test convert(Int, y) == one(Int)
        @test Float64(y) == one(Float64)
    end

    @testset "conversion array" begin
        y = Logical(1, ()->1, true)
        y.value = ones(2)
        @test convert(Vector, y) == ones(2)
    end

    @testset "indexing" begin
        y = Logical(()->1, true)
        y.value = 1
        @test y[1] == 1
    end


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
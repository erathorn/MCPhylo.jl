@testset "SubstitutionModels.jl" begin
    @testset "Restriction" begin
        target = ([-0.9191450300180579 -0.7071067811865475; 0.3939192985791677 -0.7071067811865476], 
                  [-1.0, 0.0], 
                  [-0.7615773105863908 0.7615773105863908; -0.4242640687119285 -0.9899494936611665], 
                  2.380952380952381)
        result = Restriction([0.3, 0.7], [0.0])
        @testset "Decomposition of Q matrix" begin
            @test result[1] ≈ target[1]
            @test result[2] ≈ target[2]
            @test result[3] ≈ target[3]
            @test result[4] ≈ target[4]  
            @test_throws AssertionError MCPhylo.Restriction(ones(3) / 3, [0.1])
        end

        @testset "Calculate transition matrix" begin
            U, D, Uinv, mu = result
            Q = U * diagm(D) * Uinv
            P = exp(Q * mu)
            @test P ≈ MCPhylo.calculate_transition(Restriction, 1.0, mu, 1.0, U, Uinv, D, [])
            @test [1.0 0.0; 0.0 1.0] == MCPhylo.calculate_transition(Restriction, 1.0, mu, 1.0e-12, U, Uinv, D, [])
            @test [1.0 0.0; 0.0 1.0] == MCPhylo.calculate_transition(Restriction, 1.0, mu, 1.0e-12, [1im 2im; 0im 3im], [1im 2im; 2im 3im], D, [])  
        end
    end  
    
    @testset "JC" begin
        target = ([0.7071067811865476 -0.408248290463863 -0.5773502691896257; -0.7071067811865475 -0.40824829046386313 -0.577350269189626; 0.0 0.816496580927726 -0.5773502691896258],
                  [-1.0, -0.9999999999999994, 1.1102230288049144e-16], 
                  [0.7071067811865476 -0.7071067811865475 0.0; -0.4082482904638629 -0.408248290463863 0.816496580927726; -0.5773502691896256 -0.5773502691896257 -0.5773502691896258], 
                  1.5)
        result = JC([0.15, 0.45, 0.4], [0.0])
        @testset "Decomposition of Q matrix" begin
            @test result[1] ≈ target[1]
            @test result[2] ≈ target[2]
            @test result[3] ≈ target[3]
            @test result[4] ≈ target[4]
        end  

        @testset "Calculate transition matrix" begin
            U, D, Uinv, mu = result
            Q = U * diagm(D) * Uinv
            P = exp(Q * mu)
            @test P ≈ MCPhylo.calculate_transition(JC, 1.0, mu, 1.0, U, Uinv, D, [])
            min_result = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
            @test min_result == MCPhylo.calculate_transition(JC, 1.0, mu, 1.0e-12, U, Uinv, D, [])
            @test min_result == MCPhylo.calculate_transition(JC, 1.0, mu, 1.0e-12, [0.0im 2.0im 0.0im; 0.0im 2.0im 0.0im; 0.0im 2.0im 0.0im], [1.0im 2.0im; 1.0im 2.0im], D, [])
            max_result = [0.2 0.2 0.2; 0.2 0.2 0.2; 0.2 0.2 0.2]
            @test max_result == MCPhylo.calculate_transition(JC, 1.0, mu, 101.1, [0.0im 2.0im 0.0im; 0.0im 2.0im 0.0im; 0.0im 2.0im 0.0im], [1.0im 2.0im; 1.0im 2.0im], D, [1, 1, 1, 1, 1])
        end
    end

    @testset "GTR" begin
        target = ([-0.06919383208280398 0.5545230380047614 0.5773502691896261; -0.7729814547853464 -0.6912032799070211 0.5773502691896257; 0.6306440233282372 -0.46340287673658875 0.5773502691896257],
                  [-1.3622649268468585, -0.29662234328312154, -5.551115123125783e-17], 
                  [-0.1434321592442829 -0.6409264859501984 0.7843586451944811; 0.8837782123919058 -0.44064564120807803 -0.443132571183828; 0.8660254037844387 0.34641016151377513 0.5196152422706631],
                  0.6028137161614641)
        result = GTR([0.5, 0.2, 0.3], [0.060325906174435326, 0.48940696298364417, 2.4502671308419206])
        @test result[1] ≈ target[1]
        @test result[2] ≈ target[2]
        @test result[3] ≈ target[3]
        @test result[4] ≈ target[4] 
    end

    @testset "freeK" begin
        target = ([-0.9913312346241385 -0.7071067811865475; 0.131386389167909 -0.7071067811865476], 
                  [-0.2853036708835002, -6.938893903907228e-18], [-0.8906959139221835 0.8906959139221834; 
                  -0.16549879465230727 -1.2487147677207877], 
                  3.5050372710007514)
        result = freeK([0.0], [0.03338775337571049, 0.2519159175077897, 0.8202684796606095, 2.8944278494558904])
        @test result[1] ≈ target[1]
        @test result[2] ≈ target[2]
        @test result[3] ≈ target[3]
        @test result[4] ≈ target[4]   
    end

    @testset "setmatrix" begin
        result  = MCPhylo.setmatrix([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        @test result == [0.0 0.1 0.2 0.4; 0.1 0.0 0.3 0.5; 0.2 0.3 0.0 0.6; 0.4 0.5 0.6 0.0]
    end
end
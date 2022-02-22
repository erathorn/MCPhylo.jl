@testset "Rates.jl" begin
    @testset "discrete_gamma_rates" begin
        target_vec = [0.03338775337571049, 0.2519159175077897, 0.8202684796606095, 2.8944278494558904]
        @test discrete_gamma_rates(0.5, 0.5, 4) == target_vec
        
        target_vec = [0.029077754761925846, 0.2807145371399754, 0.92477306511421, 2.7654346429838887]
        @test discrete_gamma_rates(0.5, 0.5, 4, method=:median) == target_vec
    end
    
    @testset "median_boundaries" begin
        target_vec = [0.024746651492520606, 0.23890238006213632, 0.7870290171790321, 2.353525844604548]
        @test MCPhylo.median_boundaries(0.5, 0.5, 4) == target_vec
    end

    @testset "mean_boundaries" begin
        target_vec = [0.10153104426762159, 0.45493642311957283, 1.3233036969314669, 1.3018825e-315]
        @test MCPhylo.mean_boundaries(0.5, 0.5, 4) ≈ target_vec
    end
end
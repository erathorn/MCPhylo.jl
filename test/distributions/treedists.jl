@testset "exponentialBL" begin
    exp_dist = exponentialBL(1.0)
    @test exp_dist.scale == 1.0
    @test ismissing(exp_dist.constraints)
    const_dict = Dict(:exc => [(["E", "F"], ["G"])], 
                      :mono => [["A", "B"], ["C", "D"], ["E", "F"]], 
                      :not_mono => [["C", "D", "E"]])
    exp_dist = exponentialBL(1.0, const_dict)
    @test exp_dist.scale == 1.0
    @test exp_dist.constraints == const_dict
end

@testset "CompoundDirichlet" begin
    compound_dirichlet = CompoundDirichlet(1.0,1.0,0.100,1.0)
    @test compound_dirichlet.alpha == 1.0
    @test compound_dirichlet.beta == 0.1
    @test compound_dirichlet.a == 1.0
    @test compound_dirichlet.c == 1.0
    @test ismissing(compound_dirichlet.constraints)
    const_dict = Dict(:exc => [(["E", "F"], ["G"])], 
                      :mono => [["A", "B"], ["C", "D"], ["E", "F"]], 
                      :not_mono => [["C", "D", "E"]])
    compound_dirichlet = CompoundDirichlet(1.0,1.0,0.100,1.0, const_dict)
    @test compound_dirichlet.alpha == 1.0
    @test compound_dirichlet.beta == 0.1
    @test compound_dirichlet.a == 1.0
    @test compound_dirichlet.c == 1.0
    @test compound_dirichlet.constraints == const_dict
end

@testset "TreeDistribution" begin
    ldist = exponentialBL(1.0)
    const_dict = Dict(:exc => [(["E", "F"], ["G"])], 
                      :mono => [["A", "B"], ["C", "D"], ["E", "F"]], 
                      :not_mono => [["C", "D", "E"]])
    tree_dist = TreeDistribution(ldist, const_dict)
    @test tree_dist.length_distr.scale == 1.0
    @test ismissing(tree_dist.length_distr.constraints)
    @test tree_dist.topology_distr.constraint_dict == const_dict
    @test isa(tree_dist.topology_distr, UniformConstrained)

    tree_dist = TreeDistribution(ldist)
    @test tree_dist.length_distr.scale == ldist.scale
    @test ismissing(tree_dist.length_distr.constraints)
    @test isa(tree_dist.topology_distr, UniformTopology)


    tree_dist = TreeDistribution(const_dict)
    @test isa(tree_dist.length_distr, UniformBranchLength)
    @test tree_dist.topology_distr.constraint_dict == const_dict
    @test isa(tree_dist.topology_distr, UniformConstrained)

    tree_dist = TreeDistribution()
    @test isa(tree_dist.length_distr, UniformBranchLength)
    @test isa(tree_dist.topology_distr, UniformTopology)
end
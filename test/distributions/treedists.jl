@testset "exponentialBL" begin
    exp_dist = exponentialBL(1.0)
    @test exp_dist.scale == 1.0
    @test ismissing(exp_dist.constraints)
    
    
    tree = ParseNewick(
        "(((0:0.110833,1:0.0137979)10:0.146124,(2:0.197891,(3:0.132967,(4:0.0378759,5:0.089252)11:0.101833)12:0.184301)
         13:0.0450774)14:0.335725,6:0.153197,(7:0.0216218,(8:0.0781687,9:0.120419)15:0.0209114)16:0.0209771);",
    )

    @test sum(logpdf.(Exponential(1.0), get_branchlength_vector(tree))) == logpdf(exp_dist, tree)
    @test all(-ones(17) .== gradlogpdf(exp_dist, tree)[2])
    
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


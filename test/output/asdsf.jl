@testset "asdsf.jl" begin  
    result = ASDSF("./output/trees1.nwk", "./output/trees2.nwk"; show_progress=false)
    target = 0.0802412477781256
    @test isapprox(result[end], 0.0802,atol= 0.0001)

    trees1 = newick.(ParseNewick("./output/trees1.nwk"))
    trees2 = newick.(ParseNewick("./output/trees2.nwk"))
    result = ASDSF(trees1, trees2, show_progress=false)
    @test isapprox(result[end], 0.0802,atol= 0.0001)
end
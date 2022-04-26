@testset "Felsenstein" begin

    ntax, nchar, gap, miss, symbols, df, langs =
        MCPhylo.ParseNexus("./likelihood/primates.nex")


    primates_tree = ParseNewick(
        "(((Tarsius_syrichta:0.0510942,(Lemur_catta:0.0136013,Homo_sapiens:0.0370755)12:0.0343822)
        13:0.224569,(Pan:0.0712342,Gorilla:0.03754)14:0.0295151)15:0.0768634,((Pongo:
        0.020513,Hylobates:0.159117)16:0.239429,Macaca_fuscata:0.454752)
        17:0.0902988,((M_mulatta:0.0644278,M_fascicularis:0.318016)
        18:0.015879,(M_sylvanus:0.100663,Saimiri_sciureus:0.0112774)19:0.2727)20:0.0448203);",
    )

    df = MCPhylo.datafortree(df, langs, primates_tree, symbols, gap, miss, log_space=false)
    pden = ones(4) / 4
    pd = PhyloDist(primates_tree, pden, [1.0], [1.0], JC)

    @test logpdf(pd, df) ≈ -8677.360274116634
end

@testset "Felsenstein grad" begin
    tree = ParseNewick(
        "(((0:0.110833,1:0.0137979)10:0.146124,(2:0.197891,(3:0.132967,(4:0.0378759,5:0.089252)11:0.101833)12:0.184301)
         13:0.0450774)14:0.335725,6:0.153197,(7:0.0216218,(8:0.0781687,9:0.120419)15:0.0209114)16:0.0209771);",
    )
    ntax, nchar, gap, miss, symbols, df, langs =
        MCPhylo.ParseNexus("./likelihood/simudata.nex")
    df = MCPhylo.datafortree(df, langs, tree, symbols, gap, miss, log_space=false)
    pden = ones(4) / 4
    pd = PhyloDist(tree, pden, [1.0], [1.0], JC)
    ll1, grad1 = gradlogpdf(pd, df)
    ll = -738.7363911756175
    grad = [
        -56.26844975912472, -37.95476779372888, 48.07427747703012, -13.522428720894336, -6.167080619234006, -19.397053119809236, -13.271503408059957, 22.692064516207992, -3.1448307825700947, -10.36830870394157, -6.793157740726832, -32.06664637250233, -3.747254760539362, 5.75131083892449, -2.174093469729038, 106.60905891037972, 57.1951304897995]
    
    @test ll ≈ ll1

    @test all(grad .≈ grad1)
end

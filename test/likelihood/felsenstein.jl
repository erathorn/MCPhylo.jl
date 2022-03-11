@testset "Felsenstein" begin
    
    ntax, nchar, gap, miss, symbols, df, langs = MCPhylo.ParseNexus("./likelihood/primates.nex")


    primates_tree = ParseNewick(
    "(((Tarsius_syrichta:0.0510942,(Lemur_catta:0.0136013,Homo_sapiens:0.0370755)12:0.0343822)
    13:0.224569,(Pan:0.0712342,Gorilla:0.03754)14:0.0295151)15:0.0768634,((Pongo:
    0.020513,Hylobates:0.159117)16:0.239429,Macaca_fuscata:0.454752)
    17:0.0902988,((M_mulatta:0.0644278,M_fascicularis:0.318016)
    18:0.015879,(M_sylvanus:0.100663,Saimiri_sciureus:0.0112774)19:0.2727)20:0.0448203);")

    df = MCPhylo.datafortree(df, langs, primates_tree, symbols, gap, miss)
    pden = ones(4)/4
    pd = PhyloDist(primates_tree, pden, [1.0], [1.0], JC)

    @test logpdf(pd, df) ≈ -8677.360274116634
end
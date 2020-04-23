
#=
tester:
- Julia version: 1.3.0
- Author: erathorn
- Date: 2019-05-07
=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Random
Random.seed!(42)

# support for csv
mt, df = make_tree_with_data("LangData/Dravidian.cc.phy.nex"); # load your own nexus file


po = post_order(mt);
for node in po
    node.data = df[:,:,node.num]
    node.scaler = zeros(1,size(node.data, 2))
end

my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :df => df,
  :nnodes => size(df)[3],
  :nbase => size(df)[1],
  :nsites => size(df)[2],
);



# model setup
model =  Model(
    df = Stochastic(3, (mtree, pi, rates, nnodes, nbase, nsites) -> PhyloDist(mtree, pi, rates, nbase, nsites, nnodes), false, false),
    pi = Logical( (mypi) -> mypi[1], false),
    mypi = Stochastic(1, () -> Dirichlet(2,1)),
    mtree = Stochastic(MCPhylo.Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    rates = Logical(1,(mymap, av) -> [av[convert(UInt8,i)] for i in mymap],false),
    mymap = Stochastic(1,() -> Categorical([0.25, 0.25, 0.25, 0.25]), false),
    av = Stochastic(1,() -> Dirichlet([1.0, 1.0, 1.0, 1.0]), false)
     )

# intial model values
inits = [ Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :mypi=> rand(Dirichlet(2, 1)),
    :df => my_data[:df],
    :nnodes => my_data[:nnodes],
    :nbase => my_data[:nbase],
    :nsites => my_data[:nsites],
    :mymap=>ones(3132),
    :av => [1,1,1,1]
    ),
    ]

scheme = [PNUTS(:mtree),
          #Slice(:mypi, 0.05, Univariate)
          SliceSimplex(:mypi)
          ]

setsamplers!(model, scheme);

# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 5, burnin=1,thin=1, chains=1, trees=true)

sim = mcmc(sim, 5000, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "t_dev", 5)
blv = get_branchlength_vector(mt)
using BenchmarkTools
using Zygote



"""
https://github.com/FluxML/Zygote.jl/issues/292
"""






function rtest(mt, blv, pi_, df)#::Tuple{Float64, Vector{Float64}}
    po = post_order(mt)

    nrec(y) =  MCPhylo.FelsensteinFunction(po, pi_, ones(3), df, 231, y)


    r = Zygote.pullback(nrec, blv)

    r[1],r[2](1.0)
end
pi_ = 0.8
rtest(mt, blv, pi_, df)





@Juno.profiler rtest(mt, blv, 0.8, df)



@benchmark rtest($mt, $blv, $pi_, $df)
@code_warntype MCPhylo.FelsensteinFunction(po, 0.8, ones(3), df, 231, blv)



f = open("local/newick.nwk", "r")
cont = readlines(f)
close(f)

nwkstring = cont[1]
nwkstring ="(Swedish_0:0.1034804,(Welsh_N_0:0.1422432,(Sardinian_N_0:0.02234697,(Italian_0:0.01580386,Rumanian_List_0:0.03388825):0.008238525):0.07314805):0.03669193,(((Marathi_0:0.04934081,Oriya_0:0.02689862):0.1193376,Pashto_0:0.1930713):0.05037896,Slovenian_0:0.0789572):0.03256979):0.5;"

second = "(Swedish_0:0.1034804,(Welsh_N_0:0.1422432,(Sardinian_N_0:0.02234697,(Italian_0:0.01580386,Rumanian_List_0:0.03388825)4:0.008238525)3:0.07314805)2:0.03669193)1:0.5;"
third = "((7:1.3617744475378415,(((23:1.3578093871930184,((((4:1.1897145449514694,9:0.9473887258190631)50:0.024595723457505162,((6:0.6809830402989178,17:0.8391326271724254)39:0.037847800097581395,12:1.0133687105755458)54:0.0023614313480008887)47:0.05008495299934567,25:0.8013674349404084)35:0.002179794638277786,30:1.0883567489369095)43:0.0058579708118923785)41:0.041823014378176364,((((22:1.1282669419601847,(29:1.0733609714223427,(((13:1.407796339125349,(2:0.9838809588544728,((((((((3:1.0167379679785453,(26:0.8513611941582522,8:1.0886099701304586)49:0.05642143607200655)53:0.043562588912790474,27:1.3562300261749494)57:0.058317250514426454,14:0.763418604810457)32:0.00034547344565899913,16:0.7036019944861908)38:0.03233364221332136,19:1.2334710473942765)51:0.05174746527330903,11:0.569372830853781)42:0.0122961251490966,5:1.3767993212940894)31:0.002711310068338279,21:0.8410783970100124)36:0.015656829987849504)44:0.003262971620266519)46:0.028764753346208587,15:0.6460091027873318)55:0.006922861722649104,20:0.8596200420375146)52:0.00792936416004344)48:0.020301195524972324)34:0.002727301199160799,28:0.6127844654782495)33:0.01528521687267695,1:1.192783515313331)58:0.0012483788566816625,18:1.3191757241584297)37:0.00036386004548859724)45:0.04129660447802948,10:0.7281758965361236)56:0.007723059516815107)40:0.005572416964961834,24:1.2162342198219793)59:0.5;"
mt = nwk_parser(second)
mt.children
small = "(Sardinian_N_0:0.02234697,(Italian_0:0.01580386,Rumanian_List_0:0.03388825)4:0.008238525):0.5;"
mt = second_parser([i for i in second][1:end-1])
newick(mt[1])
mt = nwk_parser(third)
newick(mt)

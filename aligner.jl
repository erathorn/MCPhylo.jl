


include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo
using Distributions
import Distributions: _logpdf, logpdf
using LinearAlgebra
using Random
using StatsFuns
Random.seed!(42)
include("asjp_dolgo_map.jl")
include("alignfunctions.jl")



fn = "LangData/ielex-northeuralex-intersection-with-asjp.tsv.asjp"
data, alphabet, lang_dict = read_cognate_data(fn)


mt = MCPhylo.create_tree_from_leaves_bin([k for k in keys(lang_dict)], 2)


l2c_dict = letter2class(mapping_classes)
i2l_dict = Dict{Int64, Char}(value=>key for (key, value) in pairs(alphabet))
rmap = get_classes(mapping_classes)
mv = maximum(values(rmap))
rmap[(-1, -1)] = mv+1
rmap[(-2, -2)] = mv+2

mymap, maparr = set_map(i2l_dict,l2c_dict, rmap)
my_data = Dict{Symbol, Any}(

  :data => data,
  :langs => size(data, 1),
  :concs => size(data, 2),
  :wsize => size(data, 3),
  :maparr => maparr,
  :nchars => 37,
  :mt => mt,
  :ts => size(mt)[1],
  :lang_dict => lang_dict
);



model =  Model(
    data = Stochastic(3,
    (emp, gx, a, r, mt) ->
		Alignment(emp, gx, gx,a,r, mt, my_data[:langs], my_data[:concs], my_data[:wsize], my_data[:lang_dict], 0.0),false),
	gx = Logical(1, (gxr) -> log.(gxr), false),
    gxr = Stochastic(1, ()->Dirichlet(37, 1.0)),
	emp = Logical(2, (emr, gxr) ->indices2values(my_data[:maparr], emr, gxr,my_data[:nchars]), false),
	emr = Stochastic(1, () -> Dirichlet(68, 1)),
	a = Stochastic( () -> Beta(1,1)),
	r = Stochastic( () ->Beta(1,1)),
	mt = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:ts]+1, true)
	    )

# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(

 :data => my_data[:data],
 :gxr => rand(Dirichlet(37,1)),
 :emr => rand(Dirichlet(68,1)),
 :a => rand(),
 :r => rand(),
 :μi => randn(),
 :σ => rand(),
 :μ => randn(),
 :mtime => rand(666),
 :mt => mt
 )
]


scheme = [SliceSimplex(:gxr),
		SliceSimplex(:emr),
		Slice(:mt, 0.05, Multivariate),
		RWM(:mt, 1),
		Slice([:a, :r], 0.05, Multivariate),
]


setsamplers!(model, scheme);

sim = mcmc(model, my_data, inits, 500, burnin=100,thin=2, chains=1, trees=true)
sim = mcmc(sim, 1000, trees=true)

to_file(sim, "alg", 2)

using Serialization
serialize("alph.jls",alphabet)
serialize("maparr.jls",maparr)
serialize("lang_dict.jls",lang_dict)
serialize("data.jls",data)
#serialize("alph.jls",alphabet)
#serialize("my_data.jls", my_data)

#write("al.jls", sim)
#my_data_2 = open(deserialize, "untracked_files/Alignment/my_data.jls", "r")
#read("al.jls", ModelChains)

#using Serialization
using CSV
#alphabet = open(deserialize, "untracked_files/Alignment/alph.jls", "r")
#log1 = CSV.read("untracked_files/Alignment/alignment_mcmc_params_1.log")
#data = open(deserialize, "untracked_files/Alignment/data.jls", "r")
#maparr = open(deserialize, "untracked_files/Alignment/maparr.jls", "r")
#f = open("untracked_files/Alignment/alignment_mcmc_mytrees_1.nwk", "r")
#nwk_raw = readlines(f)
#close(f)
#lang_dict = open(deserialize, "untracked_files/Alignment/lang_dict.jls", "r")

alphabet = open(deserialize, "alph.jls", "r")
log1 = CSV.read("algparams_1.log")
data = open(deserialize, "data.jls", "r")
maparr = open(deserialize, "maparr.jls", "r")
f = open("algmytrees_1.nwk", "r")
nwk_raw = readlines(f)
close(f)
lang_dict = open(deserialize, "lang_dict.jls", "r")



nchars = 37
langs = 30
ts = 59
wsize = 16
concs = 190
tra = log1[!,[:a, :r]]
emrnames = names(log1)[63:130]
emr = log1[!,emrnames]
gxrnames = names(log1)[132:168]
gxr = log1[!,gxrnames]


include("b_cubed.jl")



using LightGraphs
using DataFrames

gs, gs_conc, gs_langs = gs_cognate_data("LangData/ielex-northeuralex-intersection-with-asjp.tsv.asjp")



rd = eval_alignemnt(Array(tra), Array(emr), Array(gxr), maparr, nwk_raw,lang_dict,data,
                    500, 515, 1, langs, concs, wsize, q=0.98)

precs = []
recs = []
fs = []

#using Clustering

for gen in keys(rd)

	test_dict = Dict()
	for (concept, conc_ind) in gs_conc
		tempt = Dict()
		dm = rd[gen][:,:,conc_ind]


		g = DiGraph(dm)

		labels, _ = label_propagation(g)
		for (lang, lang_ind) in lang_dict
			if labels[lang_ind] in keys(tempt)
				push!(tempt[labels[lang_ind]], lang)
			else
				tempt[labels[lang_ind]]=[lang]
			end
		end
		test_dict[concept] = tempt
	end
	p, r, f = calc_f_score(gs_dict, test_dict)
	push!(precs, p)
	push!(recs, r)
	push!(fs, f)
end

mean(mean.(precs))
mean(mean.(recs))
mean(mean.(fs))




gs_dict = Dict{String,Dict{String, Array{String}}}()

for (concept, index_c) in gs_conc
	tempt_dict = Dict{String, Array{String}}()
	for (language, index_l) in gs_langs
		if gs[index_l,index_c] != ""
			if gs[index_l,index_c] in keys(tempt_dict)
				push!(tempt_dict[gs[index_l,index_c]], language)
			else
				tempt_dict[gs[index_l,index_c]]=[language]
			end
		end
	end
	gs_dict[concept] = tempt_dict
end

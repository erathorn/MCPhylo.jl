


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
		Alignment(emp, gx, gx,a,r, mt, my_data[:langs], my_data[:concs], my_data[:wsize], my_data[:lang_dict]),false),
	gx = Logical(1, (gxr) -> log.(gxr), false),
    gxr = Stochastic(1, ()->Dirichlet(37, 1.0)),
	emp = Logical(2, (emr, gxr) ->indices2values(my_data[:maparr], emr, gxr,my_data[:nchars]), false),
	emr = Stochastic(1, () -> Dirichlet(68, 1)),
	a = Stochastic( () -> Beta(1,1)),
	r = Stochastic( () -> Exponential()),
	mt = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:ts]+1, true)
	    )

# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(

 :data => my_data[:data],
 :gxr => rand(Dirichlet(37,1)),
 :emr => rand(Dirichlet(68,1)),
 :a => rand(),
 :r => rand(),

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

sim = mcmc(model, my_data, inits, 50, burnin=10,thin=2, chains=1, trees=true)

to_file(sim, "alg", 1)

using Serialization
serialize("alph.jls",alphabet)
serialize("my_data.jls", my_data)
alphabet = open(deserialize, "untracked_files/Alignment/alph.jls", "r")
write("al.jls", sim)
#my_data_2 = open(deserialize, "untracked_files/Alignment/my_data.jls", "r")
read("al.jls", ModelChains)
using CSV
log1 = CSV.read("untracked_files/Alignment/alignment_mcmc_params_1.log")
data = open(deserialize, "untracked_files/Alignment/data.jls", "r")
maparr = open(deserialize, "untracked_files/Alignment/maparr.jls", "r")
f = open("untracked_files/Alignment/alignment_mcmc_mytrees_1.nwk", "r")
nwk_raw = readlines(f)
close(f)


nchars = 37
langs = 30
ts = 59
wsize = 16
concs = 190
tr = log1[!,[:a, :r]]
emrnames = names(log1)[3:70]
emr = log1[!,emrnames]
gxrnames = names(log1)[132:168]
gxr = log1[!,gxrnames]

using DataFrames
eval_alignemnt(Array(tr), Array(emr), Array(gxr), maparr, nwk_raw, 5, 1, langs, concs, wsize)

function eval_alignemnt(tr::Array{Float64,2}, emr::Array{Float64,2}, gxr::Array{Float64,2}, maparr::Array{Int64,2},
	 nwk_raw::Vector{String}, lang_dict::Dict{String, Int64},
	 ngens::Int64, thin::Int64, langs::Int64, concs::Int64, wsize::Int64)
	outt = zeros(langs, langs, concs)
	@inbounds @views for i in 1:thin:ngens
		tri = Vector(tr[i, :])
		emri = Vector(emr[i, :])
		gxri = Vector(gxr[i, :])
		emp::Array{Float64} = indices2values(maparr, emri, gxri, nchars)
		gx::Vector{Float64} = log.(gxri)
		nwk = nwk_raw[i]
		tree = nwk_parser(strip(nwk))
		#timemat, ll = mto_distance_matrix_i(tree)
		#mappedll = [parse(Int64, node.name) for node in  ll]
		#mapper = Dict(val=>ind for (ind, val) in enumerate(mappedll))

		alg = Alignment(emp, gx, gx,tri[1],tri[2], tree , langs, concs, wsize, lang_dict)

		outt .+= cognate_statements(alg, data)
	end
	outt ./= length(1:thin:ngens)
	outt
end




tree = nwk_parser(strip(nwk_raw[50]))



function nwk_parser(newick_string)
    tree = second_parser(newick_string[1:end-1])#[1]

    MCPhylo.set_binary!(tree)
    MCPhylo.number_nodes!(tree)
    tree
end


function second_parser(nwkstring)

    ms = []
    parts = split(nwkstring, ")")

	if length(parts) == 1
	   descendants, label = [], nwkstring
   else
	   if !startswith(nwkstring, '(')
		   throw("ups")
	   end
	   descendants = parse_sibblings(join(parts[1:end-1].*")")[2:end-1])
	   label = parts[end]
	end

   name, leng = name_and_length(label)

   node = Node()
   node.name = name
   node.inc_length = leng
   for ch in descendants
	   add_child!(node,ch)
   end
   node

end

function name_and_length(ms)

    name, le = split(ms, ":")
    if length(name) == 0
        name = "no_name_"*string(le)
    end
    le = parse(Float64, le)
    return name, le
end


function parse_sibblings(nwkstring)
	bracket_level = 0
	current = []
	children = []

	for c in nwkstring * ","
		if c == ',' && bracket_level == 0
			push!(children, join(current))
			current = []
		else
			if c == '('
				bracket_level += 1
			elseif c == ')'
				bracket_level -= 1
			end
			push!(current, c)
		end
	end
	return second_parser.(children)
end

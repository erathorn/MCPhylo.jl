


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


mutable struct Alignment <: DiscreteMatrixDistribution
        emp::AbstractArray{Float64,2}
        px::AbstractArray{Float64}
        py::AbstractArray{Float64}
		a::Float64
		r::Float64
		timemat::AbstractArray{Float64}
        langs::Int64
        concs::Int64
		wsize::Int64
end

Base.minimum(d::Alignment) = -Inf
Base.maximum(d::Alignment) = +Inf
Base.size(d::Alignment) = (d.langs, d.concs, d.wsize)

function logpdf(d::Alignment, x::AbstractArray{T, 3})::Float64 where T<:Real
	calcl_ll(convert.(Int64,x), d.emp, d.px, d.py, d.a, d.r, d.timemat)
end


fn = "LangData/ielex-northeuralex-intersection-with-asjp.tsv.asjp"
data, alphabet = read_cognate_data(fn)


mt = MCPhylo.create_tree_from_leaves_bin([string(i) for i in 1:30], 2)


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
  :ts => size(mt)[1]

);



model =  Model(
    data = Stochastic(3,
    (emp, gx, a, r, timemat) ->
		Alignment(emp, gx, gx,a,r, timemat , my_data[:langs], my_data[:concs], my_data[:wsize]),false),
	gx = Logical(1, (gxr) -> log.(gxr), false),
    gxr = Stochastic(1, ()->Dirichlet(37, 1.0)),
	emp = Logical(2, (emr, gxr) ->indices2values(my_data[:maparr], emr, gxr,my_data[:nchars]), false),
	emr = Stochastic(1, () -> Dirichlet(68, 1)),
	a = Stochastic( () -> Beta(1,1)),
	r = Stochastic( () -> Beta(1,1)),
	timemat = Logical(2, (mt) -> mto_distance_matrix(mt), false),
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
		#Slice(:mtime, 0.05, Multivariate),
		#RWM(:mtime, 0.5),
		Slice(:mt, 0.05, Multivariate),
		RWM(:mt, 1),
		Slice([:a, :r], 0.05, Multivariate),
]


setsamplers!(model, scheme);

sim = mcmc(model, my_data, inits, 5, burnin=2,thin=1, chains=1, trees=true)

to_file(sim, "alg", 1)

using Serialization
serialize("alph.jls",alphabet)
serialize("my_data.jls", my_data)
alphabet = open(deserialize, "alph.jls", "r")
write("al.jls", sim)
my_data = open(deserialize, "my_data.jls", "r")
read("al.jls", ModelChains)
using CSV
log1 = CSV.read("alignmentparams_1.log")

tr = log1[!,[:a, :r]]
emrnames = names(log1)[3:70]
emr = log1[!,emrnames]
gxrnames = names(log1)[132:168]
gxr = log1[!,gxrnames]
outt = similar(my_data[:data], Float64).=0.0
for i in 1:20000
	tri = Vector(tr[i, :])
	emri = Vector(emr[i, :])
	gxri = Vector(gxr[i, :])
	emp = indices2values(my_data[:maparr], emri, gxri, my_data[:nchars])
	gx = log.(gxri)
	alg = Alignment(emp, gx, gx,tri[1],tri[2], timemat , my_data[:langs], my_data[:concs], my_data[:wsize])
	outt .+= cognate_statements(alg, data)
end
outt ./= 20000

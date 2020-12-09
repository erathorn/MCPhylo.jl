using CSV, DataFrames, DataFramesMeta

include("phoible_features.jl")

env_vars = DataFrame!(CSV.File("data/env_vars.csv"))
earthradius = 6371

sample = get_phoible_sample(phoible, 100)
sample_langs = unique(sample.Glottocode)

phoible_select = select(sample, ([:Glottocode,:LanguageName,:SegmentClass,:raisedLarynxEjective]))
env_select = select(env_vars, ([:glottocode, :latitude, :longitude, :altitude]))

glottolist_env = env2[!,:glottocode]
glottolist_phoible = ph2[!,:Glottocode]

ph3 = filter(x -> x.Glottocode in sample_langs, env_select)

consonants = @where(ph3, :SegmentClass .== "consonant")
c = unique(consonants)

ej = @where(c, :raisedLarynxEjective .== "+")

ej_list = ej[!,:Glottocode]
ej_list = Vector(ej_list)

# get all the glottocodes on that list to be represented only by + in the main table, c
for e in ej_list
   c[(c[:Glottocode] .== e),:raisedLarynxEjective]="+"
end

d = unique(c)
unique!(d, :Glottocode)
sort!(d, :Glottocode)
ejectives = @where(d, :raisedLarynxEjective .== "+")

glotto_d = Vector(d[!,:Glottocode])
env3 = filter(x -> x.glottocode in glotto_d, env2)
sort!(env3, :glottocode)

ejectives_str = convert(Array{Any}, d.raisedLarynxEjective)

# preparing the data for logistic regression with {0,1} coding

replace!(ejectives_str, "-" => 0)
replace!(ejectives_str, "+" => 1)
replace!(ejectives_str,"+,-" => 1)

d[:ejective] = ejectives_str

rename!(d,:Glottocode => "glottocode")

dc = innerjoin(d, env3, on=:glottocode)

select!(dc, Not([:SegmentClass, :raisedLarynxEjective]))
points = hcat(dc.longitude,dc.latitude)
nlang = length(dc.ejective)

logistic_data = select(dc,[:ejective,:altitude])

evec = Vector{Float64}(dc.ejective)
avec = Vector{Float64}(dc.altitude)

altitude_divided = avec ./100
dc[:altdiv] = altitude_divided

logregdf = DataFrame(y = evec, x = altitude_divided)

symmetricdf = deepcopy(dc)

replace!(symmetricdf.ejective, 0 => -1)

mpty = Matrix{Float64}(undef, nlang, nlang)
dmat = pairwise!(mpty, Haversine(earthradius), points, dims=1)

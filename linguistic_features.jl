using CSV, DataFrames, DataFramesMeta, Distances

languages = DataFrame!(CSV.File("data/languages.csv"))
l = select(languages, :ID, :Name, :Latitude, :Longitude, :Glottocode, :Genus, :ISO_codes)

geni = l[!,:Genus]
unique_genus = unique(geni)

vals = DataFrame!(CSV.File("data/values.csv"))
feature_names = DataFrame!(CSV.File("data/parameters.csv"))
v = select(vals, :Language_ID, :Parameter_ID, :Value)
fnames = select(feature_names, :ID, :Name)
rename!(v, (:Language_ID => :ID))

# all features
dc = innerjoin(v, l, on=:ID)

wo_vals = @where(v, :Parameter_ID .== "81A")
wo_languagecodes = wo_vals[!, :ID]

wo = filter(x -> x.ID in wo_languagecodes, l)

wo_data = innerjoin(wo_vals, wo, on=:ID)
dropmissing!(wo_data)

wo_values = wo_data[!,:Value]

points = hcat(wo_data.Longitude, wo_data.Latitude)
earthradius = 6357
nlang = length(wo_data.Value)
mat = Matrix{Float64}(undef, nlang, nlang)
wo_dmat = pairwise!(mat, Haversine(earthradius), points, dims=1)

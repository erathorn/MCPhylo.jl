using CSV, DataFrames, DataFramesMeta, Distances

languages = DataFrame!(CSV.File("data/languages.csv"))
languages_df = select(languages, :ID, :Name, :Latitude, :Longitude, :Glottocode, :Genus, :ISO_codes, :Samples_200)
vals = DataFrame!(CSV.File("data/values.csv"))
vals = select(vals, :Language_ID, :Parameter_ID, :Value)
rename!(vals, (:Language_ID => :ID))

# all features
dc = innerjoin(vals, languages_df, on=:ID)

"""
function get_wals_features: Args
    - a DataFrame of language data which must include the following column names:
     :ID, :Parameter_ID, :Value
    - the WALS IDs of the features you want to retrieve, as a vector of strings
"""

# I'm just doing this to fix a minor data inconvenience in the sample languages
va = @where(vals, :Parameter_ID .== "81A")
valsix = @where(va, :Value .== 6)
kxo = @where(valsix, :ID .== "kxo")

sample_languages = @where(languages_df, :Samples_200 .== true)
kxo2 = @where(languages_df, :ID .== "kxo")
kxo2 = select(kxo2, :ID, :Genus, :Latitude, :Longitude)
sample_languages = select(sample_languages, :Genus, :ID, :Latitude, :Longitude)

langs = outerjoin(sample_languages, kxo2, on=["ID", "Genus", "Latitude", "Longitude"])
sort!(langs, :ID)
dropmissing!(langs)

df_with_vals = innerjoin(langs, vals, on=:ID)
sort!(df_with_vals, :ID)

unique(df_with_vals.Genus)

"""
get_wals_features sets all missing values to -10.
"""

function get_wals_features(data::DataFrame, features::Vector{String})
    all_languages = sort!(unique(data.ID))
    nlang = length(all_languages)
    nfeatures = length(features)
    new_arr = Array{Int64,2}(undef, nfeatures, nlang)
    for (feature_idx, feature) in enumerate(features)
        features_df = @where(data, :Parameter_ID .== feature)
        languages_to_vals = Dict(features_df.ID[i] => features_df.Value[i] for i in 1:length(features_df.ID))
        for (language_idx, language) in enumerate(all_languages)
            if haskey(languages_to_vals, language)
                new_arr[feature_idx, language_idx] = languages_to_vals[language]
            else
                new_arr[feature_idx, language_idx] = -10
            end
        end
    end
    return new_arr
end

# alternative solution to missings
#data_array[ismissing.(data_array)] .= -1

feature_list = ["81A", "27A", "13A", "49A"]

data_array = get_wals_features(df_with_vals, feature_list)

function create_distance_matrix(data::DataFrame, radius::Int64)
    points = hcat(data.Longitude, data.Latitude)
    nlang = length(data.ID)
    mpty = Matrix{Float64}(undef, nlang, nlang)
    distance_matrix = pairwise!(mpty, Haversine(earthradius), points, dims=1)
    return distance_matrix
end

earthradius = 6357
sample_dmat = create_distance_matrix(langs, earthradius)

"""
array_to_dict converts the data array from get_wals_features to a dictionary.
"""

function array_to_dict(X::Array{Int64,2}, features::Array{String,1})
    nfeatures, nlang = size(X)
    d = Dict{String,Array}()
    for feature in 1:features
        setindex!(d, X[i,:], feature)
    end
    return d
end

count

function extract_nvals(X::Array{Int64,2})
    nvals = Vector{Int64}()
    nfeatures, nlang = size(X)
    for i in 1:nfeatures
        vals = count(x -> x ≠ -10, (unique(X[i,:])))
        push!(nvals, vals)
    end
    return nvals
end

extract_nvals(data_array)

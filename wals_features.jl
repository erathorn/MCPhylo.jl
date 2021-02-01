using CSV, DataFrames, DataFramesMeta, Distances

languages = DataFrame!(CSV.File("data/murawaki_data.csv"))
cols = names(languages)

function delete_pidgins_sl!(d::DataFrame)
    sl = @where(d, :genus .== "Sign Languages")
    pc = @where(d, :genus .== "Pidgins and Creoles")
    sl_and_pc = vcat(sl.wals_code, pc.wals_code)
    df = filter(x -> x ∉ sl_and_pc, d)
    return df
end

d = delete_pidgins_sl!(languages)

""" do this only if using direct-from-WALS data (not murawaki_data)

vals = DataFrame!(CSV.File("data/values.csv"))
vals = select(vals, :Language_ID, :Parameter_ID, :Value)
rename!(vals, (:Language_ID => :ID))

dc = innerjoin(vals, languages_df, on=:ID)
dc_list = unique(dc.ID)
dc_langs = filter(x -> x.ID in dc_list, languages_df)
"""
# use this for Murawaki's data

word_order = ["81A Order of Subject, Object and Verb"]

function get_feat_array(d::DataFrame, features::Vector{String})
    nlang = length(d.wals_code)
    langs = d.wals_code
    nfeat = length(features)
    arr = Array{Int64,2}(undef, nfeat, nlang)
    for (feature_idx, feature) in enumerate(features)
        val_list = d[!,feature]
        vals = unique(skipmissing(val_list))
        vals = sort!(vals)
        val_idx = Vector(1:length(vals))
        dict = Dict(vals[i] => val_idx[i] for i in 1:length(vals))
        for (lang_idx, lang) in enumerate(d.wals_code)
            if ismissing(val_list[lang_idx])
                arr[feature_idx, lang_idx] = -10
            else
                arr[feature_idx, lang_idx] = dict[val_list[lang_idx]]
            end
        end
    end
    return arr
end

data_array = get_feat_array(d, word_order)

"""
function get_wals_features: Args
    - a DataFrame of language data which must include the following column names:
     :ID, :Parameter_ID, :Value
    - the WALS IDs of the features you want to retrieve, as a vector of strings
get_wals_features sets all missing values to -10.
"""

function get_wals_features(data::DataFrame, features::Vector{String})
    all_languages = sort!(unique(data.ID))
    nlang = length(all_languages)
    nfeatures = length(features)
    new_arr = Array{Int64,2}(undef, nfeatures, nlang)
    for (feature_idx, feature) in enumerate(features)
        features_df = @where(data, :Parameter_ID .== feature)
        langs = features_df.ID
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

function create_distance_matrix(data::DataFrame, radius::Int64)
    points = hcat(data.longitude, data.latitude)
    nlang = length(data.wals_code)
    mpty = Matrix{Float64}(undef, nlang, nlang)
    distance_matrix = pairwise!(mpty, Haversine(earthradius), points, dims=1)
    return distance_matrix
end

earthradius = 6357
dmat = create_distance_matrix(d, earthradius)

function extract_nvals(X::Array{Int64,2})
    nvals = Vector{Int64}()
    nfeatures, nlang = size(X)
    for i in 1:nfeatures
        vals = count(x -> x ≠ -10, (unique(X[i,:])))
        push!(nvals, vals)
    end
    return nvals
end

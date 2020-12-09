using CSV, DataFrames, DataFramesMeta

phoible = DataFrame!(CSV.File("data/phoible.csv"))

function get_phoible_sample(data::DataFrame, sample_size::Int64)
    all_lang = sort!(unique(data.Glottocode))
    sample_lang = rand(all_lang, sample_size)
    new = filter(row -> row.Glottocode in sample_lang, data)
    return new
end

function binary_phoible_features(data::DataFrame, features::Vector{Symbol}, segmentclass::String)
    select_features = select(data, :Glottocode, :SegmentClass, features)
    new_data = @where(select_features, :SegmentClass .== segmentclass)
    all_lang = sort!(unique(data.Glottocode))
    nlang = length(all_lang)
    nfeat = length(features)
    lang_arr = Array{Union{Int64,Missing},2}(undef,nfeat,nlang)
    for (feature_idx, feature) in enumerate(features)
        pres = filter(feature => x -> x .== "+", phoible)
        lang_pres = sort!(unique(pres.Glottocode))
        abs = filter(feature => x -> x .== "-", phoible)
        lang_abs = sort!(unique(abs.Glottocode))
        for (lang_idx,lang) in enumerate(all_lang)
            if in(lang, lang_pres)
                lang_arr[feature_idx,lang_idx] = 1
            elseif in(lang, lang_abs)
                lang_arr[feature_idx,lang_idx] = -1
            else
                lang_arr[feature_idx,lang_idx] = -10
            end
        end
    end
    return lang_arr
end

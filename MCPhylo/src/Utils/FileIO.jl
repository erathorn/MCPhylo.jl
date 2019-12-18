
function to_file(model::ModelChains, outpath::AbstractString, thin)

    for run in 1:size(model.value,3)

        df = DataFrame(model.value[:,:,run])
        rename!(df, Symbol.(model.names))
        if isassigned(model.trees, 1)

            tdf = DataFrame(model.trees[:,:,run])
            to_file(df, tdf, outpath, string(run), thin)
        else
            to_file(df, outpath, string(run))
        end
    end

end
function to_file(df::DataFrame, outpath::AbstractString, run::AbstractString, thin)
    insertcols!(df,1, it=1:nrow(df))
    df[!, 1] .*= thin
    CSV.write(string(outpath, "params_"*run*".log"), df, writeheader=true, delim="\t")

end
function to_file(df::DataFrame, tdf::DataFrame, outpath::AbstractString, run::AbstractString, thin)

    insertcols!(df,1, it=1:nrow(df))
    df[!, 1] .*= thin
    CSV.write(string(outpath, "params_"*run*".log"), df, writeheader=true, delim="\t")
    io = open(string(outpath, "mytrees_"*run*".nwk"), "w")
    for x = 1:length(tdf[:,1])
        write(io, tdf[x,:][1])
        write(io, "\n")
    end
    close(io)
end

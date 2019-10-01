
function to_file(model::ModelChains, outpath::AbstractString)
    v = model.value
    df = DataFrame(model.value[:,:,1])

    names!(df, Symbol.(model.names))
    #insert!(df, 1, 1:nrow(df), :it)
    insertcols!(df,1, it=1:nrow(df))
    CSV.write(string(outpath, "params.log"), df, writeheader=true, delim="\t")
    if isassigned(model.trees, 1)
        io = open(string(outpath, "mytrees.nwk"), "w")
        for t in model.trees
            write(io, t)
            write(io, "\n")
        end
        close(io)
    end
end

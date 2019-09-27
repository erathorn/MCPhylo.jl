
function to_file(model::ModelChains, outpath::AbstractString)
    v = model.value
    df = DataFrame(model.value[:,:,1])
    names!(df, Symbol.(model.names))
    CSV.write(string(outpath, "params.csv"), df, writeheader=true, delim="\t")
    if isdefined(model.trees, 1)
        io = open(string(outpath, "mytrees.nwk"), "w")
        for t in model.trees
            write(io, t)
            write(io, "\n")
        end
        close(io)
    end
end

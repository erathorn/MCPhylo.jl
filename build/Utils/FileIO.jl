
"""
    to_file(model::ModelChains, outpath::AbstractString)

This function writes the results of the MCMC runs into files. The destination of
the files is specified using `outpath`. It will create a files for each chain. A
`params_x.log` file storing each parameter sample. In this case `x` specifies the
index of the chain. The file is compatible with MCMC analysis tools like `Tracer`
(http://tree.bio.ed.ac.uk/software/tracer/). If in addition trees are sampled,
they are stored in newick format in a file called `trees_x.nwk`, where `x` again
specifies the index of the respective chain.
"""
function to_file(model::ModelChains, outpath::AbstractString)

    for run in 1:size(model.value,3)

        df = DataFrame(model.value[:,:,run])
        rename!(df, Symbol.(model.names))
        thin = model.range.step
        if isassigned(model.trees, 1)

            tdf = DataFrame(model.trees[:,:,run])
            to_file(df, tdf, outpath, string(run), thin)
        else
            to_file(df, outpath, string(run), thin)
        end
    end

end

function to_file(df::DataFrame, outpath::AbstractString, run::AbstractString, thin::Int64)
    insertcols!(df,1, :it=>1:nrow(df))
    df[!, 1] .*= thin
    CSV.write(string(outpath, "params_"*run*".log"), df, writeheader=true, delim="\t")

end

function to_file(df::DataFrame, tdf::DataFrame, outpath::AbstractString, run::AbstractString, thin::Int64)

    insertcols!(df,1, :it=>1:nrow(df))
    df[!, 1] .*= thin
    CSV.write(string(outpath, "params_"*run*".log"), df, header=true, delim="\t")
    io = open(string(outpath, "trees_"*run*".nwk"), "w")
    for x = 1:length(tdf[:,1])
        write(io, tdf[x,:][1])
        write(io, "\n")
    end
    close(io)
end

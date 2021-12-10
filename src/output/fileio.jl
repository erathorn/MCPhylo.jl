#################### File I/O ####################
"""
    Base.read(name::AbstractString, ::Type{T}) where {T<:AbstractChains}

Read a chain from an external file.

Returns an `AbstractChains` subtype read from an external file.

* `name` : file to read or write. Recommended convention is for the file name to be specified with a `.jls` extension.

* `T` : chain type to read.

* `c` : chain to write.
"""
function Base.read(name::AbstractString, ::Type{T}) where {T<:AbstractChains}
    c = open(deserialize, name, "r")
    isa(c, T) || throw(TypeError(:open, "read(\"$name\", $T)", T, c))
    c
end
"""
    Base.write(name::AbstractString, c::AbstractChains)

Write a chain to an external file.

Returns a written external file containing a subtype.

* `name` : file to read or write. Recommended convention is for the file name to be specified with a `.jls` extension.

* `T` : chain type to read.

* `c` : chain to write.
"""
function Base.write(name::AbstractString, c::AbstractChains)
    open(file -> serialize(file, c), name, "w")
end

"""
    readcoda(output::AbstractString, index::AbstractString)

Read MCMC sampler output generated in the CODA format by OpenBUGS. The function only retains those sampler iterations at which all model parameters were monitored.

Returns a `Chains` object containing the read sampler output.

* `output` : text file containing the iteration numbers and sampled values for the model parameters.

* `index` : text file containing the names of the parameters, followed by the first and last rows in which their output can be found in the `output` file.


"""
function readcoda(output::AbstractString, index::AbstractString)
    out = readdlm(output, Any)
    ind = readdlm(index, Any)

    firstind = ind[:, 2]
    firstiter = out[firstind, 1]
    lastind = ind[:, 3]
    lastiter = out[lastind, 1]

    thin = Int((lastiter[1] - firstiter[1]) / (lastind[1] - firstind[1]))
    window = maximum(firstiter):thin:minimum(lastiter)
    startind = firstind .+ (first(window) .- firstiter) ./ step(window)
    stopind = lastind .- (lastiter .- last(window)) ./ step(window)

    names = AbstractString[ind[:, 1]...]

    value = Array{Float64}(undef, length(window), length(names))
    for i = 1:size(value, 2)
        inds = Int(startind[i]):Int(stopind[i])
        value[:, i] = out[inds, 2]
    end

    Chains(value, start = first(window), thin = step(window), names = names)
end


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
        thin = model.range.step
        if isassigned(model.trees, 1)
            to_file(model.value[:, :, run], model.trees[:,:,run], outpath, string(run), thin, model.names, model.tree_names)
        else
            to_file(model.value[:, :, run], outpath, string(run), thin, model.names)
        end
    end

end

function to_file(df::Array, outpath::String, run::String, thin::Int64, names::Array)
    io = open(string(outpath, "params_"*run*".log"),"w")
    write(io, "it\t")
    join(io, names, "\t")
    write(io, "\n")
    for (ind,row) in enumerate(eachrow(df))
        write(io, "$(ind*thin)\t")
        join(io, row, "\t")
        write(io, "\n")
    end
    close(io)
end


function to_file(df::Array, tdf::Array, outpath::AbstractString,
                run::AbstractString, thin::Int64, names::Array, tree_names::Array{A}) where A <: AbstractString

    to_file(df, outpath, run, thin, names)
    for (ind, n_name) in enumerate(tree_names)
        io = open(string(outpath, "trees_"*n_name*"_"*run*".nwk"), "w")
        for x = 1:length(tdf[:,ind])
            write(io, tdf[x,ind])
            write(io, "\n")
        end
        close(io)
    end
end

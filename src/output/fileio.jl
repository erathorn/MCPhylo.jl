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

# JSON.lower(c::AbstractChains) = [p.x, p.y]
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
  for i in 1:size(value, 2)
    inds = Int(startind[i]):Int(stopind[i])
    value[:, i] = out[inds, 2]
  end

  Chains(value, start=first(window), thin=step(window), names=names)
end

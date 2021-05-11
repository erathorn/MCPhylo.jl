#################### Model Initialization ####################
"""
    setinits!(m::Model, inits::Dict{Symbol})

Set the initial values of stochastic model nodes.

Returns the model with stochastic nodes initialized and the iter field set equal to 0.

* `m` : model with nodes to be initialized.

* `inits` : initial values for stochastic model nodes. Dictionary keys and values should be given for each stochastic node.


"""
function setinits!(m::Model, inits::Dict{Symbol})
  m.hasinputs || throw(ArgumentError("inputs must be set before inits"))
  m.iter = 0
  for key in keys(m, :dependent)
    node = m[key]
    if isa(node, AbstractStochastic) || isa(node, AbstractTreeStochastic)
      haskey(inits, key) || throw(ArgumentError("missing initial value for node : $key"))
      setinits!(node, m, inits[key])
    else
      setinits!(node, m)
    end
  end
  m.hasinits = true
  m
end
"""
    setinits!(m::Model, inits::Vector{V} where V<:Dict{Symbol})

"""
function setinits!(m::Model, inits::Vector{V} where V<:Dict{Symbol})
  n = length(inits)
  m.states = Array{ModelState}(undef, n)
  for i in n:-1:1
    setinits!(m, inits[i])
    m.states[i] = ModelState(unlist(m), deepcopy(gettune(m)))
  end
  m
end
"""
    setinputs!(m::Model, inputs::Dict{Symbol})

Set the values of input model nodes.

* `m` : model with input nodes to be assigned.

* `inputs` : values for input model nodes. Dictionary keys and values should be given for each input node.

Returns the model with values assigned to input nodes.
"""
function setinputs!(m::Model, inputs::Dict{Symbol})
  for key in keys(m, :input)
    haskey(inputs, key) ||
      throw(ArgumentError("missing inputs for node : $key"))
    isa(inputs[key], AbstractDependent) &&
      throw(ArgumentError("inputs cannot be Dependent types"))
    m.nodes[key] = deepcopy(inputs[key])
  end
  m.hasinputs = true
  m
end


"""
    function setsamplers!(m::Model, samplers::Vector{V} where V<:Sampler)

Set the samplers for the Stocastic nodes of a given Model.

Returns the Model with updated samplers.

* `m` : Model to update.

* `samplers` : block-specific sampelrs.

"""
function setsamplers!(m::Model, samplers::Vector{V} where V<:Sampler)
  m.samplers = deepcopy(samplers)
  for sampler in m.samplers
    sampler.targets = keys(m, :target, sampler.params)
  end
  m
end


"""
    SimulationParameters(; asdsf::Bool=false, freq::Int64=50, min_splits::Float64=0.1)

Construct a `SimulationParameters` object that defines the simulation parameters
for a MCMC simulation.

Returns a `SimulationParameters` type object.

* `asdsf`: controls if ASDSF should be calculated

* `freq`: at which interval trees are used for ASDSF

* `min_splits`: controls the default minimal splits threshold

"""
function SimulationParameters(; asdsf::Bool=false, freq::Int64=50,
                              min_splits::Float64=0.1)
  SimulationParameters(asdsf, freq, min_splits)
end


"""
    Base.show(io::IO, sp::SimulationParameters)

Prints the parameters of the Simulation.

* `io` : IO stream on which to print.

* `sp` : simulation parameters to print.
"""

function Base.show(io::IO, sp::SimulationParameters)
  print(io, "Object of type \"SimulationParameters\"\n\n")
  println("ASDSF = $(sp.asdsf)")
  println("ASDSF Frequency = $(sp.freq)")
  println("ASDSF minimal splits threshold = $(sp.min_splits)")
end

include("../MCPhylo.jl")
using StatsBase
using StatsPlots
using .MCPhylo
using Serialization

trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

"""
plot1 = Plots.plot(trees[1])
plot2 = Plots.plot(trees[1], treetype=:fan, msc=:blue, mc=:yellow, lc=:white,
           bg=:black, tipfont=(7, :lightgreen))
"""

data = rand(Normal(0,1), 5000)

my_data=Dict(:data=>data)

model = Model(
    data = Stochastic(1, (μ, σ) -> Normal(μ, σ), false),
       μ = Stochastic(()->Normal(),true),
       σ = Stochastic(()->Exponential(1), true)
)

inits = [Dict(:data => data,
            :μ => randn(),
            :σ => rand()),
       Dict(:data => data,
           :μ => randn(),
           :σ => rand())]

samplers = [NUTS(:μ),
           Slice(:σ, 0.1)]

setsamplers!(model, samplers)

# sim = mcmc(model, my_data, inits, 5000, burnin=500, thin=5, chains=2)
sim = mcmc(model, my_data, inits, 1000, burnin=100, thin=5, chains=2)

pv = plot2(sim, [:autocor, :mean, :density, :trace])

pv[4]
display(pv[4])

plot(pv..., layout=(4, 1), size=(1200, 800))
# savefig("test.pdf")











"""
  plot(c::AbstractChains, ptype::Vector{Symbol}=[:trace, :density];
       <keyword arguments>)::Array{Plots.Plot}

Function that takes a MCMC chain and creates various different plots (trace &
density by default).

# Arguments
- 'vars::Vector{String}=String[]' : specifies the variables of the chain that
                                    are plotted
. 'filename::String=""' : when given, the plots will be saved to a file
- 'fmt::Symbol'=:svg' : specifies the format of the output file
- 'nrow::Integer=3' / 'ncol::Integer=2' : Define layout of the plot window(s),
                                          i.e. how many plots on each page
- 'legend::Bool=false': Turn plot legend on / off
- 'args...': Plottype specific arguments, like the number of bins for the
             contourplot or if the barplots bars should be stacked or not.
            Check the specific plot functions below to use these arguments.
"""
function plot2(c::AbstractChains, ptype::Vector{Symbol}=[:trace, :density];
              vars::Vector{String}=String[], filename::String="",
              fmt::Symbol=:svg, nrow::Integer=3, ncol::Integer=2, args...
              )::Array{Plots.Plot}
  n = length(ptype)
  if !isempty(vars)
    indeces = check_vars(c.names, vars)
  else
    indeces = collect(1:length(c.names))
  end # if / else
  p = Plots.plot(c, indeces, ptype, args...)
  # draw(p, fmt=fmt, filename=filename, nrow=nrow, ncol=ncol)
  return p
end # plot


"""
  check_vars(sim_names::Vector{AbstractString},
             vars::Vector{String})::Vector{Int64}

--- INTERNAL ---
Helper function that returns a list of indeces that correspond to specific
variables. Only those variables get plotted in later steps.
"""
function check_vars(sim_names::Vector{AbstractString}, vars::Vector{String})::Vector{Int64}
    names = []
    for var in vars
        if endswith(var, r"\[[0-9]*\]")
            for sim_name in sim_names
                if sim_name == var
                    push!(names, sim_name)
                end # if
            end # for
        else
            for sim_name in sim_names
                if sim_name == var || occursin(Regex(var * "\\[[0-9]+\\]"), sim_name)
                  #(occursin(sim_name,r"\[[0-9]*\]") && sim_name[1 : findfirst("[", sim_name)[1]-1] == var))
                    push!(names, sim_name)
                end # if
            end # for
        end # if / else
    end # for
    unique!(names)
    indeces = []
    for name in names
      index = findfirst(isequal(name), sim_names)
      push!(indeces, index)
    end # for
    sort!(indeces)
    return indeces
end # check_vars


"""
    contourplot(c::AbstractChains; <keyword arguments>)::Vector{Plots.Plot}

Function that takes a MCMC chain and creates contourplots. If variables are
limited with 'vars' keyword argument, at least 2 variables have to be specified,
or no contourplot can be drawn.

# Arguments
- 'vars::Vector{String}'=String[] : specifies the variables of the chain that are plotted
. 'filename::String'="" : when given, the plots will be saved to a file
- 'fmt::Symbol'=:svg : specifies the format of the output file
- 'nrow::Integer'=3 / 'ncol::Integer'=2 : Define layout of the plot window(s),
                                          i.e. how many plots on each page
- 'legend::Bool=false': Turn plot legend on / off
"""
function contourplot(c::AbstractChains; vars::Vector{String}=String[],
                     filename::String="", fmt::Symbol=:svg, nrow::Integer=3,
                     ncol::Integer=2, legend::Bool=false, args...)::Vector{Plots.Plot}
 if !isempty(vars)
   indeces = check_vars(c.names, vars)
 else
   indeces = collect(1:length(c.names))
 end # if / else
 if length(indeces) == 1
  throw(ArgumentError("Contourplot requires at least 2 variables."))
 end # if
 p = Array{Plots.Plot}(undef, 1, length(indeces))
 showlegend = legend
 p = contourplot_int(c, indeces; legend=legend, args...)
 draw(p, fmt=fmt, filename=filename, nrow=nrow, ncol=ncol)
 return p
end # contourplot


#################### Plot Engines ####################
@recipe function f(c::AbstractChains, indeces::Vector{Int64}; ptype=[],
                   maxlag=round(Int, 10 * log10(length(c.range))),
                   position=:stack, trim=(0.025, 0.975))
  grid --> :dash
  gridalpha --> 0.5
  legend --> false
  legendtitle --> "Chain"

  arr = []
  for type in ptype

    if type == :autocor
      push!(arr, Autocor(c, indeces, maxlag))
    end # if

    if type == :bar
      push!(arr, Bar(c, indeces, position))
    end # if

    if type == :density
      push!(arr, Density(c, indeces, trim))
    end # if

    if type == :mean
      push!(arr, Mean(c, indeces))
    end # if

    if type == :trace
      push!(arr, Trace(c, indeces))
    end # if
  end # for
  return Tuple(arr)
end # recipe

struct Autocor; c; indeces; maxlag; end
struct Bar; c; indeces; position; end
struct Density; c; indeces; trim; end
struct Mean; c; indeces; end
struct Trace; c; indeces; end


@recipe function f(acor::Autocor)
  xguide --> "Lag"
  yguide --> "Autocorrelation"
  xlims --> (0, +Inf)
  layout --> (1, length(acor.indeces))

  nrows, nvars, nchains = size(acor.c.value)
  lags = 0:acor.maxlag
  ac = autocor(acor.c, lags=collect(lags))
  x = repeat(collect(lags * step(acor.c)), outer=[nchains])
  for (index, i) in enumerate(acor.indeces)
    subplot := index
    y = vec(ac.value[i,:,:])
    ac_group = repeat(acor.c.chains, inner=[length(lags)])
    for chain in acor.c.chains
      idxs = findall(==(chain), ac_group)
      @series begin
        title --> acor.c.names[i]
        seriestype := :line
        x[idxs], y[idxs]
      end # series
    end # for
  end # for
  primary := false
end # recipe


@recipe function f(bar::Bar)
  xguide --> "Value"
  yguide --> "Density"
  layout --> (1, length(bar.indeces))

  nrows, nvars, nchains = size(bar.c.value)
  for (index, i) in enumerate(bar.indeces)
    subplot := index
    S = unique(bar.c.value[:, i, :])
    n = length(S)
    x = repeat(S, 1, nchains)
    y = zeros(n, nchains)
    ymax = maximum(position == :stack ? mapslices(sum, y, dims=2) : y)
    for j in 1:nchains
      m = StatsBase.countmap(bar.c.value[:, i, j])
      for k in 1:n
        if S[k] in keys(m)
          y[k, j] = m[S[k]] / nrows
        end # if
      end # for
    end # for
    primary := false
    group := repeat(bar.c.chains, inner=[n])
    title --> bar.c.names[i]
    ylims --> (0.0, ymax)
    bar_position := position
    return StatsPlots.GroupedBar((vec(x), vec(y)))
  end # for
end # recipe


@recipe function f(dens::Density)
  xguide --> "Value"
  yguide --> "Density"
  ylims --> (0.0, +Inf)
  layout --> (1, length(dens.indeces))

  trim = (0.025, 0.975)
  nrows, nvars, nchains = size(dens.c.value)
  for (index, i) in enumerate(dens.indeces)
    subplot := index
    val = Array{Vector{Float64}}(undef, nchains)
    dens_group = []
    for j in 1:nchains
      qs = quantile(dens.c.value[:, i, j], [trim[1], trim[2]])
      mask = [qs[1] .<= dens.c.value[:, i, j] .<= qs[2]]
      val[j] = dens.c.value[mask[1], i, j]
      dens_group = vcat(dens_group, repeat([j], inner=sum(mask[1])))
    end # for
    for chain in dens.c.chains
      idxs = findall(==(chain), dens_group)
      @series begin
        title --> dens.c.names[i]
        seriestype := :density
        [val...;]
      end # series
    end # for
  end # for
  primary := false
end # recipe


@recipe function f(mean::Mean)
  xguide --> "Iteration"
  yguide --> "Mean"
  layout --> (1, length(mean.indeces))

  nrows, nvars, nchains = size(mean.c.value)
  val = cummean(mean.c.value)
  x = repeat(collect(mean.c.range), outer=[nchains])
  for (index, i) in enumerate(mean.indeces)
    subplot := index
    y = vec(val[:, i, :])
    mean_group = repeat(mean.c.chains, inner=[length(mean.c.range)])
    for chain in mean.c.chains
      idxs = findall(==(chain), mean_group)
      @series begin
        title --> mean.c.names[i]
        seriestype := :line
        x[idxs], y[idxs]
      end # series
    end # for
  end # for
  primary := false
end # recipe


@recipe function f(trace::Trace)
  xguide --> "Iteration"
  yguide --> "Value"
  layout --> (1, length(trace.indeces))
  legendtitle --> "Chain"
  widen --> false

  nrows, nvars, nchains = size(trace.c.value)
  x = repeat(collect(trace.c.range), outer=[nchains])
  for (index, i) in enumerate(trace.indeces)
    subplot := index
    y = vec(trace.c.value[:, i, :])
    trace_group = repeat(trace.c.chains, inner=[length(trace.c.range)])
    for chain in trace.c.chains
      idxs = findall(==(chain), trace_group)
      @series begin
        title --> trace.c.names[i]
        seriestype := :line
        x[idxs], y[idxs]
      end # series
    end # for
  end # for
  primary := false
end # recipe

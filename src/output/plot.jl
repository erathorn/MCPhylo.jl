#################### Posterior Plots ####################

#################### Generic Methods ####################
"""
    draw(p::Array{T}; fmt::Symbol=:svg, filename::String="", nrow::Integer=3,
         ncol::Integer=2, byrow::Bool=false, ask::Bool=true) where T<:Plots.Plot

--- INTERNAL ---
Helper function for the plot functions. Displays plots and - if wanted - saves
them to a file.
"""
function draw(p::Array{T}; fmt::Symbol=:svg, filename::String="",
              nrow::Integer=3, ncol::Integer=2, byrow::Bool=false,
              ask::Bool=true) where T<:Plots.Plot

  fmt in [:pdf, :png, :ps, :svg] ||
    throw(ArgumentError("unsupported draw format $fmt"))

  isexternalfile = length(filename) > 0
  addextension = isexternalfile && something(findfirst(isequal('.'), filename), 0) == 0

  pp = nrow * ncol               ## plots per page
  ps = length(p)                 ## number of plots
  np = ceil(Int, ps / pp)        ## number of pages
  mat = Array{Plots.Plot}(undef, pp)

  for page in 1:np
    if ask && page > 1 && !addextension
      println("Press ENTER to draw next plot")
      readline(stdin)
    end # if
    # add fileytype to filename if necessary
    if isexternalfile
      fname = filename
      if addextension
        fname = string(fname, '-', page, '.', fmt)
      end # if
    end # if

    nrem = ps - (page - 1) * pp
    for j in 1:pp
      if j <= nrem
        # store all plots for the page in an array
        mat[j] = p[(page -1) * pp + j]

      else
        # invisible empty plot
        mat[j] = Plots.plot(showaxis=:false, grid=:false, axis=nothing, border=nothing)
      end # if else
    end # for
    result = byrow ? permutedims(reshape(mat, ncol, nrow), [2, 1]) :
                     reshape(mat, nrow, ncol)
    # draw plot and save to file
    plots = Plots.plot(result..., layout=(nrow,ncol), widen=false,
               guidefontsize=8, titlefontsize=10)
    !isexternalfile && display(plots)
    isexternalfile && savefig(fname)
  end # for
end # draw


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
function plot(c::AbstractChains, ptype::Vector{Symbol}=[:trace, :density];
              vars::Vector{String}=String[], filename::String="",
              fmt::Symbol=:svg, nrow::Integer=3, ncol::Integer=2, args...
              )::Array{Plots.Plot}
  n = length(ptype)
  if !isempty(vars)
    indeces = check_vars(c.names, vars)
  else
    indeces = collect(1:length(c.names))
  end # if / else
  p = Vector{Plots.Plot}(undef, length(ptype))
  for i in 1:n
    # showlegend = legend && i == n
    p[i] = Plots.plot(c, indeces, type=ptype[i], args...)
  end # for
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
@recipe function f(c::AbstractChains, indeces::Vector{Int64}; type=:autocor,
                   maxlag=round(Int, 10 * log10(length(c.range))),
                   position=:stack)
  grid --> :dash
  gridalpha --> 0.5
  legend --> false
  legendtitle --> "Chain"

  if type == :autocor
    Autocor(c, indeces, maxlag)
  end # if

  if type == :bar
    Bar(c, indeces, position)
  end # if

  if type == :mean
    Mean(c, indeces)
  end # if
end # recipe


struct Autocor; c; indeces; maxlag; end
struct Bar; c; indeces; position; end
struct Density; c; indeces; trim; end
struct Mean; c; indeces; end
struct Trace; c; indeces; end


@recipe function f(autoc::Autocor)
  xguide --> "Lag"
  yguide --> "Autocorrelation"
  xlims --> (0, +Inf)
  layout --> length(autoc.indeces)

  nrows, nvars, nchains = size(autoc.c.value)
  lags = 0:autoc.maxlag
  ac = autocor(autoc.c, lags=collect(lags))
  x = repeat(collect(lags * step(autoc.c)), outer=[nchains])
  for (index, i) in enumerate(autoc.indeces)
    subplot := index
    y = vec(ac.value[i,:,:])
    ac_group = repeat(autoc.c.chains, inner=[length(lags)])
    for chain in autoc.c.chains
      idxs = findall(==(chain), ac_group)
      @series begin
        title --> autoc.c.names[i]
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
  layout --> length(bar.indeces)

  nrows, nvars, nchains = size(bar.c.value)
  for (index, i) in enumerate(bar.indeces)
    subplot := index
    S = unique(bar.c.value[:, i, :])
    n = length(S)
    x = repeat(S, 1, nchains)
    y = zeros(n, nchains)
    ymax = maximum(position == :stack ? mapslices(sum, y, dims=2) : y)
    for j in 1:nchains
      m = countmap(c.value[:, i, j])
      for k in 1:n
        if S[k] in keys(m)
          y[k, j] = m[S[k]] / nrows
        end # if
      end # for
    end # for
    group := repeat(bar.c.chains, inner=[n])
    title --> bar.c.names[i]
    ylims --> (0.0, ymax)
    bar_position := position
    StatsPlots.GroupedBar((vec(x), vec(y)))
  end # for
  primary := false
end # recipe

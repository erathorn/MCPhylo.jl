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
              vars::Vector{String}=String[], fuse::Bool=false, f_layout=nothing,
              fsize::Tuple{Number, Number}=(0, 0), args...
              )::Union{Vector{Plots.Plot}, Tuple{Vector{Plots.Plot}, Plots.Plot}}
  if !isempty(vars)
    indeces = check_vars(c.names, vars)
  else
    indeces = collect(1:length(c.names))
  end # if / else
  ilength = length(indeces)
  if :contour in ptype && ilength == 1
    filter!(e -> e ≠ :contour, ptype)
    if isempty(ptype)
      throw(ArgumentError("Drawing a contourplot requires at least 2 variables."))
    else
      @warn "Drawing a contourplot requires at least 2 variables. Only drawing other plots"
    end # if/else
  end # if
  n = length(ptype)
  p = Array{Plots.Plot}(undef, length(ptype))
  for i in 1:n
    p[i] = Plots.plot(c, indeces; ptype=ptype[i], size=(ilength * 500, 300), args...)
    if !fuse
      if n != 1
        println("Press ENTER to draw next plot")
        readline(stdin)
      end # if
      display(p[i])
    end # if
  end # for
  if fuse
    isnothing(f_layout) && (f_layout = (length(ptype), 1))
    fsize == (0, 0) && (fsize = (ilength * 500, length(ptype) * 300))
    allplots = Plots.plot(p..., layout=f_layout, size=fsize)
    display(allplots)
    return (p, allplots)
  end # if
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


#################### Plot Engines ####################
@recipe function f(c::AbstractChains, indeces::Vector{Int64}; ptype=:autocor,
                   maxlag=round(Int, 10 * log10(length(c.range))),
                   position=:stack, bins=100, trim=(0.025, 0.975))
  grid --> :dash
  gridalpha --> 0.5
  legend --> false
  legendtitle --> "Chain"
  margin --> 7mm

  arr = []
  ptype == :autocor ? push!(arr, Autocor(c, indeces, maxlag)) :
  ptype == :bar     ? push!(arr, Bar(c, indeces, position)) :
  ptype == :contour ? push!(arr, Contour(c, indeces, bins)) :
  ptype == :density ? push!(arr, Density(c, indeces, trim)) :
  ptype == :mean    ? push!(arr, Mean(c, indeces)) :
  ptype == :trace   ? push!(arr, Trace(c, indeces)) :
    throw(ArgumentError("unsupported plot type $ptype"))
  return Tuple(arr)
end # recipe

struct Autocor; c; indeces; maxlag; end
struct Bar; c; indeces; position; end
struct Contour; c; indeces; bins; end
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
    for (index, chain) in enumerate(acor.c.chains)
      idxs = findall(==(chain), ac_group)
      @series begin
        label --> string(index)
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


@recipe function f(cont::Contour)
  layout --> (1, sum(collect(1:length(cont.indeces) - 1)))

  subplot_index = 0
  nrows, nvars, nchains = size(cont.c.value)
  offset = 1e4 * eps()
  n = nrows * nchains
  for (index, i) in enumerate(cont.indeces[1:end-1])
    X = cont.c.value[:, i, :]
    qx = range(minimum(X) - offset, stop=maximum(X) + offset, length=cont.bins + 1)
    mx = map(k -> mean([qx[k], qx[k + 1]]), 1:cont.bins)
    idx = Int[findfirst(k -> qx[k] <= x < qx[k + 1], 1:cont.bins) for x in X]
    for j in cont.indeces[index+1:end]
      subplot_index += 1
      Y = cont.c.value[:, j, :]
      qy = range(minimum(Y) - offset, stop=maximum(Y) + offset, length=cont.bins + 1)
      my = map(k -> mean([qy[k], qy[k + 1]]), 1:cont.bins)
      idy = Int[findfirst(k -> qy[k] <= y < qy[k + 1], 1:cont.bins) for y in Y]
      density = zeros(cont.bins, cont.bins)
      for k in 1:n
        density[idx[k], idy[k]] += 1.0 / n
      end
      @series begin
        subplot := subplot_index
        xguide --> cont.c.names[i],
        yguide --> cont.c.names[j]
        colorbar_title --> "Density"
        title --> cont.c.names[i]
        seriestype := :contour
        mx, my, density
      end # series
    end # for
  end # for
  primary := false
end # recipe


@recipe function f(dens::Density)
  xguide --> "Value"
  yguide --> "Density"
  ylims --> (0.0, +Inf)
  layout --> (1, length(dens.indeces))

  nrows, nvars, nchains = size(dens.c.value)
  for (index, i) in enumerate(dens.indeces)
    subplot := index
    val = Array{Vector{Float64}}(undef, nchains)
    dens_group = []
    for j in 1:nchains
      qs = quantile(dens.c.value[:, i, j], [dens.trim[1], dens.trim[2]])
      mask = [qs[1] .<= dens.c.value[:, i, j] .<= qs[2]]
      val[j] = dens.c.value[mask[1], i, j]
      dens_group = vcat(dens_group, repeat([j], inner=sum(mask[1])))
    end # for
    for (index, chain) in enumerate(dens.c.chains)
      idxs = findall(==(chain), dens_group)
      @series begin
        label --> string(index)
        title --> dens.c.names[i]
        seriestype := :density
        [val...;][idxs]
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
    for (index, chain) in enumerate(mean.c.chains)
      idxs = findall(==(chain), mean_group)
      @series begin
        label --> string(index)
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
    for (index, chain) in enumerate(trace.c.chains)
      idxs = findall(==(chain), trace_group)
      @series begin
        label --> string(index)
        title --> trace.c.names[i]
        seriestype := :line
        x[idxs], y[idxs]
      end # series
    end # for
  end # for
  primary := false
end # recipe

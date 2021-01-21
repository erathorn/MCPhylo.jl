#################### Posterior Plots ####################

"""

plot(sim) ==> draw(plot(sim))
plot(sim, varnames) ==> plot selected variables
draw ==> nur für output
"""




#################### Generic Methods ####################

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
        mat[j] = Plots.plot(showaxis=:false, grid=:false)
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
end # function


function plot(c::AbstractChains, ptype::Vector{Symbol}=[:trace, :density];
              indeces::Vector{Integer}=[], filename::String="",
              fmt::Symbol=:svg, nrow::Integer=3, ncol::Integer=2,
              legend::Bool=false, args...)::Array{Plots.Plot}
  n = length(ptype)
  if !isempty(var_names)
    indeces = check_var_names(c.names, var_names)
  else
    indedes = [collect(1:length(c.names))]
  end # if / else
  p = Array{Plots.Plot}(undef, n, length(var_names)
  for i in 1:n
    showlegend = legend && i == n
    p[i, :] = plot(c, ptype[i], indeces; legend=legend, args...)
  end # for
  draw(p, fmt=fmt, filename=filename, nrow=nrow, ncol=ncol)
  return p
end # function


function check_var_names(sim_names::Vector{String}, var_names::Vector{String})::Vector{Integer}
    names = []
    for var_name in var_names
        if endswith(var_name, r"\[[0-9]*\]")
            for sim_name in sim_names
                if sim_name == var_name
                    push!(names, sim_name)
                end # if
            end # for
        else
            for sim_name in sim_names
                if sim_name == var_name || occursin(Regex(var_name * "\\[[0-9]+\\]"), sim_name)
                  #(occursin(sim_name,r"\[[0-9]*\]") && sim_name[1 : findfirst("[", sim_name)[1]-1] == var_name))
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
end # function


function plot(c::AbstractChains, ptype::Symbol, indeces::Vector{Integer}; legend::Bool=false, args...)
  ptype == :autocor      ? autocorplot(c, indeces; legend=legend, args...) :
  ptype == :bar          ? barplot(c, indeces; legend=legend, args...) :
  ptype == :contour      ? contourplot(c, indeces; args...) :
  ptype == :density      ? densityplot(c, indeces; legend=legend, args...) :
  ptype == :mean         ? meanplot(c, indeces; legend=legend, args...) :
  ptype == :mixeddensity ? mixeddensityplot(c, indeces; legend=legend, args...) :
  ptype == :trace        ? traceplot(c, indeces; legend=legend, args...) :
    throw(ArgumentError("unsupported plot type $ptype"))
end # function


#################### Plot Engines ####################

function autocorplot(c::AbstractChains, indeces::Vector{Integer};
                     maxlag::Integer=round(Int, 10 * log10(length(c.range))),
                     legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plots.Plot}(undef, length(indeces))
  pos = legend ? :right : :none
  lags = 0:maxlag
  ac = autocor(c, lags=collect(lags))
  j = 1
  for i in indeces
    # new plot creation block, based on Plots with a GR backend
    plots[j] = Plots.plot(repeat(collect(lags * step(c)), outer=[nchains]),
                          vec(ac.value[i,:,:]), seriestype=:line,
                          group=repeat(c.chains, inner=[length(lags)]),
                          xlabel="Lag", ylabel="Autocorrelation",
                          title=c.names[i], legendtitle="Chain", legend=pos,
                          xlims=(0, +Inf), grid=:dash; gridalpha=0.5)
    j += 1
  end # for
  return plots
end # function


function barplot(c::AbstractChains, indeces::Vector{Integer};
                 legend::Bool=false, position::Symbol=:stack, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plots.Plot}(undef, length(indeces))
  pos = legend ? :right : :none
  z = 1
  for i in indeces
    S = unique(c.value[:, i, :])
    n = length(S)
    x = repeat(S, 1, nchains)
    y = zeros(n, nchains)
    for j in 1:nchains
      m = countmap(c.value[:, i, j])
      for k in 1:n
        if S[k] in keys(m)
          y[k, j] = m[S[k]] / nrows
        end # if
      end # for
    end # for
    ymax = maximum(position == :stack ? mapslices(sum, y, dims=2) : y)
    # new plot creation block, based on StatsPlots with a GR backend
    plots[z] = StatsPlots.groupedbar(vec(x), vec(y), bar_position=position,
                                     group=repeat(c.chains, inner=[n]),
                                     xlabel="Value", ylabel="Density",
                                     title=c.names[i], legendtitle="Chain",
                                     legend=pos, ylims=(0.0, ymax),
                                     grid=:dash, gridalpha=0.5)
    z += 1
  end # for
  return plots
end # function


function contourplot(c::AbstractChains, indeces::Vector{Integer};
                     bins::Integer=100, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Plots.Plot[]
  offset = 1e4 * eps()        [1, [4, 5]]
  n = nrows * nchains
  for index, i in enumerate(indeces[1:end-1])
    X = c.value[:, i, :]
    qx = range(minimum(X) - offset, stop=maximum(X) + offset, length=bins + 1)
    mx = map(k -> mean([qx[k], qx[k + 1]]), 1:bins)
    idx = Int[findfirst(k -> qx[k] <= x < qx[k + 1], 1:bins) for x in X]
    for j in indeces[index+1:end]
      Y = c.value[:, j, :]
      qy = range(minimum(Y) - offset, stop=maximum(Y) + offset, length=bins + 1)
      my = map(k -> mean([qy[k], qy[k + 1]]), 1:bins)
      idy = Int[findfirst(k -> qy[k] <= y < qy[k + 1], 1:bins) for y in Y]
      density = zeros(bins, bins)
      for k in 1:n
        density[idx[k], idy[k]] += 1.0 / n
      end
      # new plot creation block, based on Plots with a GR backend
      p = Plots.plot(mx, my, density, seriestype=:contour,
                     colorbar_title="Density", xlabel=c.names[i],
                     ylabel=c.names[j])
      push!(plots, p)
    end # for
  end # for
  return plots
end  ##  function


function densityplot(c::AbstractChains, indeces::Vector{Integer};
                     legend::Bool=false,
                     trim::Tuple{Real, Real}=(0.025, 0.975), na...)
  nrows, nvars, nchains = size(c.value)
  # new list initialization
  plots = Array{Plots.Plot}(undef, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    val = Array{Vector{Float64}}(undef, nchains)
    grouping = []
    for j in 1:nchains
      qs = quantile(c.value[:, i, j], [trim[1], trim[2]])
      # keep all values between qs[1] and qs[2]
      mask = [qs[1] .<= c.value[:, i, j] .<= qs[2]]
      val[j] = c.value[mask[1], i, j]
      grouping = vcat(grouping, repeat([j], inner=sum(mask[1])))
    end # for
    # new plot creation block, based on StatsPlots with a GR backend
    plots[i] = StatsPlots.plot([val...;], seriestype=:density,
                               group = grouping, xlabel="Value",
                               ylabel="Density", title=c.names[i],
                               legendtitle="Chain", legend=pos,
                               ylims=(0.0, +Inf), grid=:dash, gridalpha=0.5)
  end # for
  return plots
end # function


function meanplot(c::AbstractChains, indeces::Vector{Integer}; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  # new list initialization
  plots = Array{Plots.Plot}(undef, nvars)
  pos = legend ? :right : :none
  val = cummean(c.value)
  for i in 1:nvars
    # new plot creation block, based on Plots with a GR backend
    plots[i] = Plots.plot(repeat(collect(c.range), outer=[nchains]),
                          vec(val[:, i, :]), seriestype=:line,
                          group=repeat(c.chains, inner=[length(c.range)]),
                          xlabel="Iteration", ylabel="Mean", title=c.names[i],
                          legendtitle="Chain", legend=pos,
                          grid=:dash, gridalpha=0.5)
  end # for
  return plots
end # function


function mixeddensityplot(c::AbstractChains, indeces::Vector{Integer};
                          barbounds::Tuple{Real, Real}=(0, Inf), args...)
  plots = Array{Plots.Plot}(undef, size(c, 2))
  discrete = MCPhylo.indiscretesupport(c, barbounds)
  barplots = plot(c, :bar; args...)
  densityplots = plot(c, :density; args...)
  for i in 1:length(discrete)
    if discrete[i] == true
      plots[i] = barplots[i]
    else
      plots[i] = densityplots[i]
    end # if / else
  end # for
  return plots
end # function


function traceplot(c::AbstractChains, indeces::Vector{Integer}; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  # new list initialization
  plots = Any[]
  pos = legend ? :right : :none
  for i in 1:nvars
    push!(plots, Plots.plot(repeat(collect(c.range), outer=[nchains]),
                            vec(c.value[:, i, :]), seriestype=:line,
                            group=repeat(c.chains, inner=[length(c.range)]),
                            xlabel="Iteration", ylabel="Value",
                            title=c.names[i], legendtitle="Chain", legend=pos,
                            widen=false, grid=:dash, gridalpha=0.5))
  end  # for
  return plots
end # function

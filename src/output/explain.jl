include("../MCPhylo.jl")
using .MCPhylo

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

sim = mcmc(model, my_data, inits, 5000, burnin=500, thin=5, chains=2)

# default "inner" layout puts plots in a row
pv = plot2(sim, [:mean])
# "inner" layout can be manipulated, but usually size has to be adjusted as well
pv = plot2(sim, [:mean], layout=(3, 1), size=(800,1500))
# throws an error, as it should for contour (when only one variable is selected)
pv = plot2(sim, [:contour], vars=["likelihood"])
# gives a warning for contourplot but shows the other ptypes
pv = plot2(sim, [:contour, :density, :mean], vars=["likelihood"], fuse=true)
# specific plot variables are passed successfully
pv = plot2(sim, [:autocor, :contour, :density, :mean, :trace],
           maxlag=10, bins=20, trim=(0.1, 0.9), legend=true)
# demonstrate the customizable "outer" layout
pv = plot2(sim, [:autocor, :bar, :contour, :density, :mean, :trace],
           fuse=true, fLayout=(2,2), fsize=(2750, 2500), linecolor=:match)
# barplot works
pv = plot2(sim, [:bar], linecolor=:match, legend=:true)
# use savefig to save as file; no draw function needed
savefig("test.pdf")










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
  p = Array{Plots.Plot}(undef, n)
  for i in 1:n
    if ptype[i] == :bar
      p[i] = bar_int(c, indeces; args...)
    else
      p[i] = plot(c, indeces; ptype=ptype[i], size=(ilength * 500, 300), args...)
    end # if / else
    if !fuse
      if n != 1
        println("Press ENTER to draw next plot")
        readline(stdin)
      end # if
      display(p[i])
    end # if
  end # for
  if fuse
    isnothing(f_layout) && (f_layout = (n, 1))
    fsize == (0, 0) && (fsize = (ilength * 500, n * 300))
    allplots = plot(p..., layout=f_layout, size=fsize)
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
                   maxlag=round(Int, 10 * log10(length(c.range))), bins=100,
                   trim=(0.025, 0.975))
  grid --> :dash
  gridalpha --> 0.5
  legend --> false
  legendtitle --> "Chain"
  legendtitlefonthalign := :left
  margin --> 7mm

  arr = []
  ptype == :autocor ? push!(arr, Autocor(c, indeces, maxlag)) :
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


function bar_int(c::AbstractChains, indeces::Vector{Int64};
                 position::Symbol=:stack, args...)
  nrows, nvars, nchains = size(c.value)
  ilength = length(indeces)
  bar_plots = Array{Plots.Plot}(undef, ilength)
  for (index, i) in enumerate(indeces)
    S = unique(c.value[:, i, :])
    n = length(S)
    x = repeat(S, 1, nchains)
    y = zeros(n, nchains)
    for j in 1:nchains
      m = StatsBase.countmap(c.value[:, i, j])
      for k in 1:n
        if S[k] in keys(m)
          y[k, j] = m[S[k]] / nrows
        end # if
      end # for
    end # for
    ymax = maximum(position == :stack ? mapslices(sum, y, dims=2) : y)
    bar_plots[index] = groupbar(vec(x), vec(y);
                                group=repeat(c.chains, inner=[n]),
                                title=c.names[i], ylims=(0.0, ymax), args...)
  end # for
  p = plot(bar_plots..., layout=(1, ilength), size=(ilength * 500, 300))
  return p
end # function


@recipe function f(cont::Contour)
  layout --> (1, sum(collect(1:length(cont.indeces) - 1)))
  legend := false
  colorbar --> :best
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
        xguide --> cont.c.names[i]
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


@userplot GroupBar

recipetype(::Val{:groupbar}, args...) = GroupBar(args)

Plots.group_as_matrix(g::GroupBar) = true

grouped_xy(x::AbstractVector, y::AbstractArray) = x, y
grouped_xy(y::AbstractArray) = 1:size(y,1), y

@recipe function f(g::GroupBar; spacing = 0)
    xguide --> "Value"
    yguide --> "Density"
    grid --> :dash
    gridalpha --> 0.5
    legend --> false
    legendtitle --> "Chain"
    legendtitlefonthalign := :left
    margin --> 7mm
    x, y = grouped_xy(g.args...)

    nr, nc = size(y)
    isstack = pop!(plotattributes, :bar_position, :dodge) == :stack
    isylog = pop!(plotattributes, :yscale, :identity) ∈ (:log10, :log)
    # the_ylims = pop!(plotattributes, :ylims, (-Inf, Inf))

    # extract xnums and set default bar width.
    # might need to set xticks as well
    xnums = if eltype(x) <: Number
        xdiff = length(x) > 1 ? mean(diff(x)) : 1
        bar_width --> 0.8 * xdiff
        x
    else
        bar_width --> 0.8
        ux = unique(x)
        xnums = (1:length(ux)) .- 0.5
        xticks --> (xnums, ux)
        xnums
    end
    @assert length(xnums) == nr

    # compute the x centers.  for dodge, make a matrix for each column
    x = if isstack
        x
    else
        bws = plotattributes[:bar_width] / nc
        bar_width := bws * clamp(1 - spacing, 0, 1)
        xmat = zeros(nr,nc)
        for r=1:nr
            bw = Plots._cycle(bws, r)
            farleft = xnums[r] - 0.5 * (bw * nc)
            for c=1:nc
                xmat[r,c] = farleft + 0.5bw + (c-1)*bw
            end
        end
        xmat
    end

    fill_bottom = if isylog
        if isfinite(the_ylims[1])
            min(minimum(y) / 100, the_ylims[1])
        else
            minimum(y) / 100
        end
    else
        0
    end
    # compute fillrange
    y, fr = isstack ? groupedbar_fillrange(y) : (y, get(plotattributes, :fillrange, [fill_bottom]))
    if isylog
        replace!( fr, 0 => fill_bottom )
    end
    fillrange := fr

    seriestype := :bar
    x, y
end

function groupedbar_fillrange(y)
    nr, nc = size(y)
    # bar series fills from y[nr, nc] to fr[nr, nc], y .>= fr
    fr = zeros(nr, nc)
    y = copy(y)
    y[.!isfinite.(y)] .= 0
    for r = 1:nr
        y_neg = 0
        # upper & lower bounds for positive bar
        y_pos = sum([e for e in y[r, :] if e > 0])
        # division subtract towards 0
        for c = 1:nc
            el = y[r, c]
            if el >= 0
                y[r, c] = y_pos
                y_pos -= el
                fr[r, c] = y_pos
            else
                fr[r, c] = y_neg
                y_neg += el
                y[r, c] = y_neg
            end
        end
    end
    y, fr
end

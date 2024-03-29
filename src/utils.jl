#################### Model Expression Operators ####################

function modelfx(literalargs::Vector{Tuple{Symbol,DataType}}, f::Function)
    modelfxsrc(literalargs, f)[1]
end

function modelfxsrc(literalargs::Vector{Tuple{Symbol,DataType}}, f::Function)
    args = Expr(:tuple, map(arg -> Expr(:(::), arg[1], arg[2]), literalargs)...)
    expr, src = modelexprsrc(f, literalargs)
    fx = Core.eval(Main, Expr(:function, args, expr))
    (fx, src)
end


function modelexprsrc(f::Function, literalargs::Vector{Tuple{Symbol,DataType}})

    m = first(methods(f).ms) # returns a mehtod datatype
    argnames = Vector{Any}(undef, m.nargs)
    argnames = Base.method_argnames(m)
    fkeys = Symbol[argnames[2:end]...]
    
    ftypes = DataType[m.sig.parameters[2:end]...]
    n = length(fkeys)

    literalinds = Int[]
    for (key, T) in literalargs
        i = findfirst(fkey -> fkey == key, fkeys)
        if i !== nothing && ftypes[i] == T
            push!(literalinds, i)
        end
    end
    nodeinds = setdiff(1:n, literalinds)

    all(T -> T == Any, ftypes[nodeinds]) ||
        throw(ArgumentError("model node arguments are not all of type Any"))

    modelargs = Array{Any}(undef, n)
    for i in nodeinds
        modelargs[i] = Expr(:call, MCPhylo.mod_value, Expr(:ref, :model, QuoteNode(fkeys[i])))
    end
    for i in literalinds
        modelargs[i] = fkeys[i]
    end
    expr = Expr(:block, Expr(:(=), :f, f), Expr(:call, :f, modelargs...))
    (expr, fkeys[nodeinds])
end

mod_value(d::AbstractVariate) = d.value
mod_value(d) = d

#################### Mathematical Operators ####################

# isprobvec(p::AbstractVector) = isprobvec(convert(Vector{Float64}, p))
isprobvec(s::AbstractVariate) = isprobvec(s.value)

cummean(x::AbstractArray) = mapslices(cummean, x, dims = 1)

function cummean(x::AbstractVector{T}) where {T<:Real}
    y = similar(x, T)
    xs = zero(T)
    @inbounds @fastmath for i = 1:length(x)
        xs += x[i]
        y[i] = xs / i
    end
    y
end


## Csorgo S and Faraway JJ. The exact and asymptotic distributions of the
## Cramer-von Mises statistic. Journal of the Royal Statistical Society,
## Series B, 58: 221-234, 1996.
function pcramer(q::Real)
    p = 0.0
    for k = 0:3
        c1 = 4.0 * k + 1.0
        c2 = c1^2 / (16.0 * q)
        p += gamma(k + 0.5) / factorial(k) * sqrt(c1) * exp(-c2) * besselk(0.25, c2)
    end
    p / (pi^1.5 * sqrt(q))
end


#################### Auxiliary Functions ####################
## pmap2 is a partial work-around for the pmap issue in julia 0.4.0 of worker
## node errors being blocked.  In single-processor mode, pmap2 calls map
## instead to avoid the error handling issue.  In multi-processor mode, pmap is
## called and will apply its error processing.

function pmap2(f::F, lsts::Vector) where F<:Function
    if (nprocs() > 1) && (length(lsts) > 1)
        pmap(f, lsts)
    else
        map(f, lsts)
    end
end

## Custom Macro to get embedded errror messages
macro async_showerr(ex)
    quote
        t = @async try
            eval($(esc(ex)))
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
        end
    end
end

"""
  assign_mcmc_work(f::Function, lsts::AbstractArray
                   )::Tuple{Vector{Tuple{Chains, Model, ModelState}},
                      Array{Float64, 2}, Vector{AbstractString}}

--- INTERNAL ---
Starts the MCMC chain generation (on multiple workers if possible) and also
starts parallel ASDSF - if possible and requested by the user.
"""
function assign_mcmc_work(
    f::Function,
    lsts::Vector,
    sp::SimulationParameters,
    conv_storage::Union{Nothing,ConvergenceStorage}
 )::Tuple{
     Vector{Tuple{Chains,Model,ModelState}},
     Array{Float64,2},
     Vector{AbstractString},
     Union{Nothing,ConvergenceStorage},
 }

    ASDSF::Bool = sp.asdsf
    statnames::Vector{String} = []
    # count the number of trees per step per chain
    tree_dim::Int64 = 0
    for i in lsts[1][1].nodes
        if isa(i[2], Stochastic{<:GeneralNode})
            push!(statnames, string("asdsf_", string(i[1])))
            tree_dim += 1
        end # if
    end # for
    
    nchains::Int64 = length(lsts)
    ntrees::Int64 = ASDSF ? floor((last(lsts[1][3]) - lsts[1][4]) / sp.freq) : 0
    results = Dict{Int64,Tuple{Chains,Model,ModelState}}()
    if ASDSF
        # set up remote channels for communication across workers for ASDSF-on-the-fly
        if nworkers() >= nchains + 1
            r_channels =
                [RemoteChannel(() -> Channel{Vector{AbstractString}}(100)) for x = 1:nchains]
        else
            r_channels =
                [RemoteChannel(() -> Channel{Vector{AbstractString}}(ntrees)) for x = 1:nchains]
        end # if / else
        ASDSF_vals::Vector{Vector{Float64}} = []
        # add necessary ASDSF related arguments to lsts for mcmc_worker function
        for (ind, lst) in enumerate(lsts)
            append!(lst, [sp.freq, r_channels[ind]])
        end # for
        # add a list to lsts array, that contains args for calculate_convergence
        push!(lsts, [sp, conv_storage, r_channels, ntrees, 1:tree_dim])
    end # if

    # set up a ProgressMeter for each chain...
    meters = [Progress(lsts[1][3][end] - lsts[1][3][1] + 1; desc="Chain $c: ", enabled=sp.verbose, offset=c-1, showspeed=true) for c in 1:nchains]
    # ... & one RemoteChannel that will be used to communicate the updates of the progress bar across parallel functions
    channel = RemoteChannel(() -> Channel{Integer}(1))
    # add the RemoteChannels to the lists that are later passed to the mcmc_worker
    for c in 1:nchains
        insert!(lsts[c], 8, (channel, c))
    end # for
    results_vec = []
    @sync begin
        # check RemoteChannel for new entries and updates the ProgressMeters
        finished_chains = 0
        @async while finished_chains < nchains
            chain = take!(channel)
            chain > 0 ? ProgressMeter.next!(meters[chain]) : finished_chains += 1
        end # while
        results_vec = pmap2(f, lsts)
    end # @sync
    ASDSF && close.(r_channels)
    if ASDSF
        asdsf_results = results_vec[end]
        append!(ASDSF_vals, asdsf_results[1])
        conv_storage = asdsf_results[2]
        stats = Array{Float64,2}(undef, length(ASDSF_vals[1]), length(ASDSF_vals))
        for i = 1:length(ASDSF_vals)
            stats[:, i] = ASDSF_vals[i]
        end # for
    else
        stats = zeros(Float64, 0, tree_dim)
    end # if/else
    println("\n")
    return [results_vec[i] for i = 1:nchains], stats, statnames, conv_storage
end # assign_mcmc_work


"""
  mcmc_or_convergence(args::AbstractArray)
    ::Union{Tuple{Chains, Model, ModelState}},
            Tuple{Vector{Vector{Float64}}, ConvergenceStorage}

--- INTERNAL ---
Used in pmap call in assign_mcmc_work to correctly dispatch arguments to either
calculate convergence statistics or mcmc chains.
"""
function mcmc_or_convergence(args::AbstractArray
                            )::Union{
                            Tuple{Chains, Model, ModelState},
                            Tuple{Vector{Vector{Float64}}, ConvergenceStorage}
                            }
    # depending on the elements of the args array either start calculate_convergence or...
    if isa(args[1], SimulationParameters)
        calculate_convergence(args...)
    # ...mcmc_worker
    else
        # if no convergence statistics are wanted, then the mcmc_worker will be called with exactly 8 elements in args...
        if length(args) == 8
            mcmc_worker!(args)
        # ...but with on-the-fly convergence statistics, we have 2 more elements (the step range & the remote channel), 
        # that we need to dispatch correctly to the mcmc_worker
        else
            mcmc_worker!(args[1:end-2], args[end-1:end]...)
        end # if/else
    end # if/else
end


ind2sub(dims, ind) = Tuple(CartesianIndices(dims)[ind])

showall(v) = showall(stdout, v)

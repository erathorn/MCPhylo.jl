"""
    parse_and_number(treestring::S)::Node where S<:AbstractString

--- INTERNAL ---
Parse a newick string and then call set_binary and number_nodes! on the
resulting tree
"""
function parse_and_number(treestring::S)::Node where S<:AbstractString
    p_tree2 = parsing_newick_string(string(treestring))
    set_binary!(p_tree2)
    number_nodes!(p_tree2)
    p_tree2
end # parse_and_number


"""
    nodnums2names(bps::Vector{Tuple}, tree::Node)
        ::Vector{Tuple{Set{String}, Set{String}}}

--- INTERNAL ---
"""
function nodnums2names(bps::Vector{Tuple}, tree::Node
                      )::Vector{Tuple{Set{String}, Set{String}}}

    treebps = Vector{Tuple{Set{String}, Set{String}}}()
    for bipart in bps
        bipartition = []
        for part in bipart
            myset = Set{String}()
            for num in part
                no = find_num(tree, num)
                if no.nchild == 0
                    push!(myset, no.name)
                end # if
            end # for
            push!(bipartition, myset)
        end # for
        push!(treebps, Tuple(bipartition))
    end # for
    treebps
end # nodnums2names


"""
    ASDSF(args::String...; freq::Int64=1, check_leaves::Bool=true,
          min_splits::Float64=0.1, show_progress::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for two or more
files containing newick representations of trees. Default frequency is 1 and
by default only trees with the same leafsets are supported. The default minimal
splits threshold is 0.1. The progress bar is activated by default.
"""
function ASDSF(args::String...; freq::Int64=1, check_leaves::Bool=true,
               min_splits::Float64=0.1, show_progress::Bool=true
               )::Vector{Float64}

    length(args) < 2 && throw(ArgumentError("At least two input files are needed."))
    splitsQueue = [Accumulator{Tuple{Set{String}, Set{String}}, Int64}()]
    splitsQueues = [Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()]
    for arg in args
        push!(splitsQueues[1], Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    iter = zip([eachline(arg) for arg in args]...)
    ASDF_vals = [zeros(Int(countlines(args[1]) / freq))]
    ASDSF_int(splitsQueue, splitsQueues, iter, 1, ASDF_vals, freq, check_leaves,
              min_splits, show_progress, basic=true)[1]
end # ASDSF


"""
    ASDSF(args::Vector{String}...; freq::Int64=1, check_leaves::Bool=true,
          min_splits::Float64=0.1, show_progress::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for two or more
Vectors containing newick representations of trees. Default frequency is 1 and
by default only trees with the same leafsets are supported. The default minimal
splits threshold is 0.1. The progress bar is activated by default.
"""
function ASDSF(args::Vector{String}...; freq::Int64=1, check_leaves::Bool=true,
               min_splits::Float64=0.1, show_progress::Bool=true
               )::Vector{Float64}

    length(args) < 2 && throw(ArgumentError("At least two input arrays are needed."))
    splitsQueue = [Accumulator{Tuple{Set{String}, Set{String}}, Int64}()]
    splitsQueues = [Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()]
    for arg in args
        push!(splitsQueues[1], Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    iter = zip(args...)
    ASDF_vals = [zeros(Int(length(iter) / freq))]
    ASDSF_int(splitsQueue, splitsQueues, iter, 1, ASDF_vals, freq, check_leaves,
              min_splits, show_progress; basic=true)[1]
end # ASDSF


"""
    ASDSF(model::ModelChains; freq::Int64=1, check_leaves::Bool=true,
          min_splits::Float64=0.1, show_progress::Bool=true
          )::Vector{Vector{Float64}}

Calculate the average standard deviation of split frequencies for the trees in
different chains in a ModelChains object. Default frequency is 1 and by default
only trees with the same leafsets are supported. The default minimal splits
threshold is 0.1. The progress bar is activated by default.
"""
function ASDSF(model::ModelChains; freq::Int64=1, check_leaves::Bool=true,
               min_splits::Float64=0.1, show_progress::Bool=true
               )::Vector{Vector{Float64}}
    tree_dims::UnitRange{Int64} = 1:size(model.trees, 2)
    splitsQueue = [Accumulator{Tuple{Set{String}, Set{String}}, Int64}() for x in tree_dims]
    splitsQueues = [Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}() for x in tree_dims]
    nchains = size(model.trees, 3)
    for i in 1:nchains
        for j in tree_dims
            push!(splitsQueues[j], Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
        end # for
    end # for
    if length(tree_dims) > 1
        trees = Array{Vector{AbstractString}, 2}(undef, size(model.trees, 1), nchains)
        for i in 1:size(model.trees, 1)
            for j in 1:nchains
                trees[i, j] = model.trees[i,:,j]
            end # for
        end # for
    end # if
    iter = zip([trees[:,c] for c in 1:nchains]...)
    ASDF_vals::Vector{Vector{Float64}} = [zeros(Int(floor(length(iter) / freq))) for x in tree_dims]
    ASDSF_int(splitsQueue, splitsQueues, iter, tree_dims, ASDF_vals, freq,
              check_leaves, min_splits, show_progress)
end # ASDSF


"""
    ASDSF(r_channels::Vector{RemoteChannel}, n_trees::Int64,
          tree_dims::UnitRange{Int64}, min_splits::Float64
          )::Vector{Vector{Float64}}

--- INTERNAL ---
Calculates - on-the-fly - the average standard deviation of split frequencies
for the trees generated by MCMC draws from a model. Takes a vector of remote
channels (where the generated trees are stored during the mcmc simulation) and
the total number of trees in each chain as arguments. The default minimal splits
threshold is 0.1.
"""
function ASDSF(r_channels::Vector{RemoteChannel{Channel{Array{AbstractString,1}}}},
               n_trees::Int64, tree_dims::UnitRange{Int64}, min_splits::Float64
               )::Vector{Vector{Float64}}

    iter = 1:n_trees
    ASDF_vals::Vector{Vector{Float64}} = [zeros(Int(n_trees)) for x in tree_dims]
    splitsQueue = [Accumulator{Tuple{Set{String}, Set{String}}, Int64}() for x in tree_dims]
    splitsQueues = [Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}() for x in tree_dims]
    nchains = length(r_channels)
    for i in 1:nchains
        for j in tree_dims
            push!(splitsQueues[j], Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
        end # for
    end # for
    ASDSF_int(splitsQueue, splitsQueues, iter, tree_dims, ASDF_vals, 1, false,
              min_splits, false; r_channels=r_channels)
end # ASDSF


"""
    ASDSF_int(splitsQueue, splitsQueues, iter, tree_dims, ASDF_vals, freq,
              check_leaves, min_splits, show_progress; r_channels=nothing
              )::Vector{Float64}

--- INTERNAL ---
Handles the computation of the Average Standard Deviation of Split Frequencies.
"""
function ASDSF_int(splitsQueue, splitsQueues, iter, tree_dims, ASDF_vals, freq,
                   check_leaves, min_splits, show_progress; r_channels=nothing,
                   basic=false)::Vector{Vector{Float64}}

    all_keys = [Set{Tuple{Set{String},Set{String}}}() for x in tree_dims]
    if show_progress
        p = Progress(length(ASDF_vals))
    end # if
    run::Int64 = 1
    outer::Float64 = 0.0
    inner::Float64 = 0.0
    for (i, line) in enumerate(iter)
        if mod(i, freq) == 0
            if !isnothing(r_channels)
                line = [take!(rc) for rc in r_channels]
            end
            for td in tree_dims
                trees = basic ? [parse_and_number(tree) for tree in line] :
                                [parse_and_number(tree[td]) for tree in line]
                check_leaves && check_leafsets(trees)
                bps = [get_bipartitions(tree) for tree in trees]
                named_bps = [nodnums2names(bp, tree) for (bp, tree) in zip(bps, trees)]
                cmds = [countmap(named_bp) for named_bp in named_bps]
                outer = 0.0
                new_splits = union([keys(cmds[t]) for t in 1:length(trees)]...)
                all_keys[td] = union(all_keys[td], new_splits)
                for split in new_splits
                    inc!(splitsQueue[td], split)
                end # for
                for split in all_keys[td]
                    inner = 0.0
                    for (t, tree) in enumerate(trees)
                        if haskey(cmds[t], split)
                            inner += 1 / length(trees)
                            inc!(splitsQueues[td][t], split)
                        end # if
                    end # for
                    if any([x[split] / run > min_splits for x in splitsQueues[td]])
                        fre = splitsQueue[td][split] / run
                        inner += (inner - fre) ^ 2
                        outer += sqrt(inner)
                    end # if
                end # for
                ASDF_vals[td][run] = outer / length(splitsQueue[td])
            end # for
            run += 1
            show_progress && ProgressMeter.next!(p)
        end # if
    end # for
    ASDF_vals
end # ASDSF_int

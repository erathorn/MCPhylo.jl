"""
    parse_and_number(treestring::S)::FNode where S<:AbstractString

--- INTERNAL ---
Parse a newick string and then call set_binary and number_nodes! on the
resulting tree
"""
function parse_and_number(treestring::S)::FNode where S<:AbstractString
    p_tree2 = parsing_newick_string(string(treestring))
    set_binary!(p_tree2)
    number_nodes!(p_tree2)
    p_tree2
end # parse_and_number


"""
    nodnums2names(bps::Vector{Tuple}, tree::GeneralNode)
        ::Vector{Tuple{Set{String}, Set{String}}}

--- INTERNAL ---
"""
function nodnums2names(bps::Vector{Tuple}, tree::GeneralNode
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
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    iter = zip([eachline(arg) for arg in args]...)
    ASDF_vals = zeros(Int(countlines(args[1]) / freq))
    asdsf_int(splitsQueue, splitsQueues, iter, ASDF_vals, freq, check_leaves,
              min_splits, show_progress)
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
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    iter = zip(args...)
    ASDF_vals = zeros(Int(length(iter) / freq))
    asdsf_int(splitsQueue, splitsQueues, iter, ASDF_vals, freq, check_leaves,
              min_splits, show_progress)
end # ASDSF


"""
    ASDSF(model::ModelChains; freq::Int64=1, check_leaves::Bool=true,
          min_splits::Float64=0.1, show_progress::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for the trees in
different chains in a ModelChains object. Default frequency is 1 and by default
only trees with the same leafsets are supported. The default minimal splits
threshold is 0.1. The progress bar is activated by default.
"""
function ASDSF(model::ModelChains; freq::Int64=1, check_leaves::Bool=true,
               min_splits::Float64=0.1, show_progress::Bool=true
               )::Vector{Float64}
               
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    l = size(model.trees, 3)
    for i in 1:l
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    iter = zip([model.trees[:,:,i] for i in 1:l]...)
    ASDF_vals = zeros(Int(length(iter) / freq))
    asdsf_int(splitsQueue, splitsQueues, iter, ASDF_vals, freq, check_leaves,
              min_splits, show_progress)
end # ASDSF


"""
    asdsf_int(splitsQueue, splitsQueues, iter, ASDF_vals, freq, check_leaves
              min_splits)::Vector{Float64}

--- INTERNAL ---
Handles the computation of the Average Standard Deviation of Split Frequencies.
"""
function asdsf_int(splitsQueue, splitsQueues, iter, ASDF_vals, freq,
                   check_leaves, min_splits, show_progress)::Vector{Float64}

    all_keys = Set{Tuple{Set{String},Set{String}}}()
    if show_progress
        p = Progress(length(ASDF_vals))
    end # if
    run::Int64 = 1
    outer::Float64 = 0.0
    inner::Float64 = 0.0
    for (i, lines) in enumerate(iter)
        if mod(i, freq) == 0
            trees = [parse_and_number(tree_string) for tree_string in lines]
            check_leaves && check_leafsets(trees)
            bps = [get_bipartitions(tree) for tree in trees]
            named_bps = [nodnums2names(bp, tree) for (bp, tree) in zip(bps, trees)]
            cmds = [countmap(named_bp) for named_bp in named_bps]
            outer = 0.0
            new_splits = union([keys(cmds[t]) for t in 1:length(trees)]...)
            all_keys = union(all_keys, new_splits)
            for split in new_splits
                inc!(splitsQueue, split)
            end # for
            for split in all_keys
                inner = 0.0
                for (t, tree) in enumerate(trees)
                    if haskey(cmds[t], split)
                        inner += 1 / length(trees)
                        inc!(splitsQueues[t], split)
                    end # if
                end # for
                if any([x[split] / run > min_splits for x in splitsQueues])
                    fre = splitsQueue[split] / run
                    inner += (inner - fre) ^ 2
                    outer += sqrt(inner)
                end # if
            end # for
            ASDF_vals[run] = outer / length(splitsQueue)
            run += 1
            show_progress && ProgressMeter.next!(p)
        end # if
    end # for
    ASDF_vals
end # asdsf_int

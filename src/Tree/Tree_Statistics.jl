function parse_and_number(treestring::S)::Node where S<:AbstractString
    p_tree2 = parsing_newick_string(string(treestring))
    set_binary!(p_tree2)
    number_nodes!(p_tree2)
    p_tree2
end # parse_and_number


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
    ASDSF(args...; freq::Int64=1, check_leaves::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for two or more
files containing newick representations of trees.
"""
function ASDSF(args::String...; freq::Int64=1, check_leaves::Bool=true
              )::Vector{Float64}
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    length(args) < 2 && throw(ArgumentError("At least two input files are needed."))
    iterator = zip([eachline(arg) for arg in args]...)
    asdsf_int(splitsQueue, splitsQueues, iterator, freq, check_leaves)
end # ASDSF


function ASDSF(args::Vector{String}...; freq::Int64=1, check_leaves::Bool=true
              )::Vector{Float64}
    length(args) < 2 && throw(ArgumentError("At least two input arrays are needed."))
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    ASDF_vals = zeros(Int(length(args[1]) / freq))
    iterator = zip(args...)
    asdsf_int(splitsQueue, splitsQueues, iterator, freq, check_leaves)
end # ASDSF


function ASDSF(args::ModelChains...; freq::Int64=1, check_leaves::Bool=true
              )::Vector{Float64}
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    length(args) < 2 && throw(ArgumentError("At least two input arrays are needed."))
    iterator = zip(args...)
    asdsf_int(splitsQueue, splitsQueues, iterator, freq, check_leaves)
end # ASDSF


function asdsf_int(splitsQueue, splitsQueues, iterator, freq, check_leaves
                  )::Vector{Float64}
    all_keys = Set{Tuple{Set{String},Set{String}}}()
    ASDF_vals = zeros(Int(length(iterator) / freq))
    p = Progress(length(ASDF_vals))
    run = 1
    for (i, lines) in enumerate(iterator)
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
                if any([x[split] / run > 0.1 for x in splitsQueues])
                    fre = splitsQueue[split] / run
                    inner += (inner - fre) ^ 2
                    outer += sqrt(inner)
                end # if
            end # for
            ASDF_vals[run] = outer / length(splitsQueue)
            run += 1
            ProgressMeter.next!(p)
        end # if
    end # for
    ASDF_vals
end # asdsf_int

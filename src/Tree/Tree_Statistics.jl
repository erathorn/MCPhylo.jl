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
    ASDSF(args...; freq::Int64=1, check_leafsets::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for two files
containing newick representations of trees.
"""
function ASDSF(args...; freq::Int64=1, check_leaves::Bool=true
              )::Vector{Float64}

    function asdsf_int(i, lines)
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
    end # asdsf_int

    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end # for
    all_keys = Set{Tuple{Set{String},Set{String}}}()
    run = 1
    if all(x-> typeof(x) == Vector{String}, args)
        ASDF_vals = zeros(Int(length(args[1]) / freq))
        args_joined = Vector{String}()
        for arg in args
            push!(args_joined, join(arg, "\n"))
        end # for
        p = Progress(length(ASDF_vals))
        for (i, lines) in enumerate(zip([eachline(IOBuffer(arg)) for arg in args_joined]...))
            asdsf_int(i, lines)
        end # for
    elseif all(x-> typeof(x) == String, args)
        ASDF_vals = zeros(Int(countlines(args[1]) / freq))
        p = Progress(length(ASDF_vals))
        for (i, lines) in enumerate(zip([eachline(arg) for arg in args]...))
            asdsf_int(i, lines)
        end # for
    else
        throw(ArgumentError("Input needs to be filenames as strings, or arrays
                            of newick strings"))
    end # if/else
    ASDF_vals



end # ASDSF

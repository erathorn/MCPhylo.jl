
"""
calculate_convergence(sm::SimulationParameters, args...)::Vector{Vector{Float64}}

Function that checks which
"""
function calculate_convergence(sm::SimulationParameters,
                           conv_storage::Union{Nothing,ConvergenceStorage},
                           args...
                          )::Tuple{Vector{Vector{Float64}},ConvergenceStorage}
    if sm.asdsf
        if isnothing(conv_storage)
            ASDSF_vals, conv_storage = ASDSF(args..., sm.min_splits)
        else
            ASDSF_vals, conv_storage = ASDSF(args..., sm.min_splits, cs=conv_storage)
        end # if/else
    end
    
    return ASDSF_vals, conv_storage
end


function PNUTS_monitor(sim::AbstractChains)
    NNIs = []
    tree_depths = []
    nchains = length(sim.model.states)
    for ch in 1:nchains
        for (ind, val) in enumerate(sim.model.states[ch].value)
            if isa(val, GeneralNode)
                push!(NNIs, sim.model.states[ch].tune[ind].moves)
                push!(tree_depths, sim.model.states[ch].tune[ind].tree_depth_trace)
            end
        end
    end

    p=Plots.plot(1:sim.model.iter, NNIs[1], title="NNI moves")
    for ch in 2:nchains
        Plots.plot!(1:sim.model.iter, NNIs[ch])
    end
    display(p)
    println("Press ENTER to draw next plot")
    readline(stdin)
    p=Plots.plot(1:sim.model.iter, tree_depths[1], title="Search depth")
    for ch in 2:nchains
        Plots.plot!(1:sim.model.iter, tree_depths[ch])
    end
    display(p)
    println("Summary Stats")
    println("Per Chain mean number of NNI moves per generation")
    mn = mean.(NNIs)
    println("Chain 1: $(mn[1])")
    println("Chain 2: $(mn[2])")
    println("Per Chain mean search tree depths per generation")
    mn = mean.(tree_depths)
    println("Chain 1: $(mn[1])")
    println("Chain 2: $(mn[2])")
    nothing
end

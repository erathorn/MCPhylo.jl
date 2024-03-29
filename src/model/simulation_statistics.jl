
"""
calculate_convergence(sm::SimulationParameters, args...)::Vector{Vector{Float64}}

Function that checks if asdsf needs to be calculated and if a ConvergenceStorage struct 
needs to be dispatched.
"""
function calculate_convergence(
    sm::SimulationParameters,
    conv_storage::Union{Nothing,ConvergenceStorage},
    args...,
)::Tuple{Vector{Vector{Float64}},ConvergenceStorage}
    if sm.asdsf
        if isnothing(conv_storage)
            ASDSF_vals, conv_storage = ASDSF(args..., sm.min_splits)
        else
            ASDSF_vals, conv_storage = ASDSF(args..., sm.min_splits, cs = conv_storage)
        end # if/else
    end

    return ASDSF_vals, conv_storage
end


function PNUTS_monitor(sim::AbstractChains)
    NNIs = []
    tree_depths = []
    accs = []
    nchains = length(sim.model.states)
    for ch = 1:nchains
        for (ind, val) in enumerate(sim.model.states[ch].value)
            if isa(val, GeneralNode)
                push!(NNIs, sim.model.states[ch].tune[ind].moves)
                push!(tree_depths, sim.model.states[ch].tune[ind].tree_depth_trace)
                push!(accs, sim.model.states[ch].tune[ind].acc_p_r)
            end
        end
    end

    p = Plots.plot(1:sim.model.iter, NNIs[1], title = "NNI moves")
    for ch = 2:nchains
        Plots.plot!(1:sim.model.iter, NNIs[ch])
    end
    display(p)
    println("Press ENTER to draw next plot")
    readline(stdin)
    p = Plots.plot(1:sim.model.iter, tree_depths[1], title = "Search depth")
    for ch = 2:nchains
        Plots.plot!(1:sim.model.iter, tree_depths[ch])
    end
    display(p)
    println("Press ENTER to draw next plot")
    readline(stdin)
    p = Plots.plot(1:sim.model.iter, accs[1], title = "acc_p_r")
    for ch = 2:nchains
        Plots.plot!(1:sim.model.iter, accs[ch])
    end
    display(p)

    println("Per Chain mean number of NNI moves per generation")
    mn = mean.(NNIs)
    for i = 1:nchains
        println("Chain $i: $(mn[i])")
    end
    println("Per Chain mean search tree depths per generation")
    mn = mean.(tree_depths)
    for i = 1:nchains
        println("Chain $i: $(mn[i])")
    end
    println("Per Chain mean acceptances per generation")
    mn = mean.(accs)
    for i = 1:nchains
        println("Chain $i: $(mn[i])")
    end
    nothing
end

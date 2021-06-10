"""
    generate_constraints(; mono::Vector{Vector{String}}=Vector{String}[],
                         not_mono::Vector{Vector{String}}=Vector{String}[],
                         exc::Vector{Vector{String}}=Vector{String}[]
                         )::Dict{Symbol, Vector{Vector{String}}}

Generate a dictionary of constraints based on the given arguments.
Mono for all minophyletic groups. not_mono for leafs that are not allowed
to form a monophletic group. And exc for a partial constraint, where one or more
leafs are not allowed to be part of a specific clade. The excluded leafs are
separated from the clade they are not a part of via a colon entry in the array.

e.g.: generate_constraints(mono=[["a", "b"], ["a", "c"]],
                           exc=["b", "c", ":", "d"])
"""
function generate_constraints(; mono::Vector{Vector{String}}=Vector{String}[],
                              not_mono::Vector{Vector{String}}=Vector{String}[],
                              exc::Vector{Vector{String}}=Vector{String}[]
                              )::Dict{Symbol, Vector{Vector{String}}}

    constraints_dict = Dict{Symbol, Vector{Vector{String}}}()
    constrainttypes = [:mono, :not_mono, :exc]
    constraints_dict[:mono] = unique(mono)
    constraints_dict[:not_mono] = unique(not_mono)
    constraints_dict[:exc] = unique(exc)
    return constraints_dict
end # generate_constraints


"""
    generate_constraints!(constraints::Dict;
                        mono::Vector{Vector{String}}=Vector{String}[],
                        not_mono::Vector{Vector{String}}=Vector{String}[],
                        exc::Vector{Vector{String}}=Vector{String}[]
                        )::Dict{Symbol, Vector{Vector{String}}}

Function that adds further constraints to an existing dictionary of constraints.
See basic generate_constraints function for more info.
"""
function generate_constraints!(constraints::Dict;
                               mono::Vector{Vector{String}}=Vector{String}[],
                               not_mono::Vector{Vector{String}}=Vector{String}[],
                               exc::Vector{Vector{String}}=Vector{String}[]
                               )::Dict{Symbol, Vector{Vector{String}}}

    constraints2 = generate_constraints(mono=mono, not_mono=not_mono, exc=exc)
    # make sure entries are unique for each constraint category
    for key in keys(constraints)
        constraints[key] = union(constraints[key], constraints2[key])
    end # for
    return constraints
end # generate_constraints


"""
    generate_constraints(filename::String)

Function that creates a dictionary of constraints, based on a txt file with a
specific format (each line one constraint), i.e.:

mono        ,a,b,c
mono        ,d,e
not_mono    ,c,d
exc         ,a,b,:,d

See basic generate_constraints function for more info.
"""
function generate_constraints(filename::String)
    mono, not_mono, exc = parse_constraints(filename)
    generate_constraints(mono=mono, not_mono=not_mono, exc=exc)
end # generate_constraints


"""
    generate_constraints!(constraints::Dict{Symbol, Vector{Vector{String}}},
                          filename::String)

Function that adds further constraints to an existing dictionary of constraints,
based on a txt file with a specific format (each line one constraint). See basic
generate_constraints function on a file for more info.
"""
function generate_constraints!(constraints::Dict{Symbol, Vector{Vector{String}}}, filename::String)
    mono, not_mono, exc = parse_constraints(filename)
    generate_constraints!(constraints; mono=mono, not_mono=not_mono, exc=exc)
end # generate_constraints!


"""
    parse_constraints(filename::String)

--- INTERNAL ---
Helper function that reads constraints out of a file
"""
function parse_constraints(filename::String)
    mono, not_mono, exc = [Vector{String}[] for _ = 1:3]
    for line in eachline(filename)
        splitted_line = split(line, ",")
        splitted_line[1] = strip(splitted_line[1])
        splitted_line[1] == "mono" ? push!(mono, splitted_line[2:end]) :
        splitted_line[1] == "not_mono" ? push!(not_mono, splitted_line[2:end]) :
        splitted_line[1] == "exc" ? push!(exc, splitted_line[2:end]) :
        @warn "Unsupported constraint type $(splitted_line[1]). Allowed types are 'mono', 'not_mono' and 'exc'"
    end # for
    return (mono, not_mono, exc)
end # parse_constraints


function topological(tree::N, constraints::Dict) where N<:GeneralNode
    for key in keys(constraints)
        lca = find_lca(tree, constraints[key])
        lca.root && return false
    end
    true
end


# topological constraints fallback
function topological(tree::N, constraints::Missing) where N<:GeneralNode
    true
end

"""
    generate_constraints(
        ; mono::Vector{Vector{S}}=Vector{String}[],
        not_mono::Vector{Vector{S}}=Vector{String}[],
        exc::Vector{Tuple{Vector{S}, Vector{S}}}=Vector{Tuple{Vector{String}, Vector{String}}}[]
    )::Dict{Symbol, Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}}
        where S<:AbstractString

Generate a dictionary of constraints based on the given arguments.
Mono for all minophyletic groups. not_mono for leafs that are not allowed
to form a monophletic group. And exc for a partial constraint, where one or more
leafs are not allowed to be part of a specific clade.

e.g.: generate_constraints(mono=[["a", "b"], ["a", "c"]],
                           exc=[(["b", "c"],["d"]])
"""
function generate_constraints(
    ; mono::Vector{Vector{S}}=Vector{String}[],
    not_mono::Vector{Vector{S}}=Vector{String}[],
    exc::Vector{Tuple{Vector{S}, Vector{S}}}=Vector{Tuple{Vector{String}, Vector{String}}}[]
)::Dict{
    Symbol, Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}} where S<:AbstractString

    constraints_dict = Dict{Symbol, Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}}()
    constrainttypes = [:mono, :not_mono, :exc]
    constraints_dict[:mono] = unique(mono)
    constraints_dict[:not_mono] = unique(not_mono)
    constraints_dict[:exc] = unique(exc)
    return constraints_dict
end # generate_constraints


"""
    generate_constraints!(
        constraints::Dict;
        mono::Vector{Vector{S}}=Vector{String}[],
        not_mono::Vector{Vector{S}}=Vector{String}[],
        exc::Vector{Tuple{Vector{S}, Vector{S}}}=Vector{Tuple{Vector{String}, Vector{String}}}[]
    )::Dict{Symbol, Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}} where S<:AbstractString

Function that adds further constraints to an existing dictionary of constraints.
See basic generate_constraints function for more info.
"""
function generate_constraints!(constraints::Dict;
                               mono::Vector{Vector{S}}=Vector{String}[],
                               not_mono::Vector{Vector{S}}=Vector{String}[],
                               exc::Vector{Tuple{Vector{S}, Vector{S}}}=Vector{Tuple{Vector{String}, Vector{String}}}[]
                               )::Dict{Symbol, Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}} where S<:AbstractString


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
exc         ,a,b:d,e

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
function generate_constraints!(constraints::Dict, filename::String)
    mono, not_mono, exc = parse_constraints(filename)
    generate_constraints!(constraints; mono=mono, not_mono=not_mono, exc=exc)
end # generate_constraints!


"""
    parse_constraints(filename::String)

--- INTERNAL ---
Helper function that reads constraints out of a file
"""
function parse_constraints(filename::String
                          )::Tuple{Vector{Vector{S}}, Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}} where S<:AbstractString

    exc = Vector{Tuple{Vector{AbstractString}, Vector{AbstractString}}}()
    mono, not_mono = [Vector{String}[] for _ = 1:2]
    for line in eachline(filename)
        line = split(line,":")
        line = [split(x, ",") for x in line]
        line[1][1] = strip(line[1][1])
        line[1][1] == "mono" ? push!(mono, line[1][2:end]) :
        line[1][1] == "not_mono" ? push!(not_mono, line[1][2:end]) :
        line[1][1] == "exc" ? push!(exc, (line[1][2:end], line[2])) :
        @warn "Unsupported constraint type $(line[1][1]). Allowed types are 'mono', 'not_mono' and 'exc'"
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

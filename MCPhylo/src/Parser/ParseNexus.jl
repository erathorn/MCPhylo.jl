

"""
    ParseNexus(filename::String)

This function parses a NEXUS file which stores the input for the MCMC compuation.
The file should follow the conventions used for MrBayes.
"""
function ParseNexus(filename::String)
    open(filename, "r") do file
        global content = readlines(file)
    end # do

    if content[1] != "#NEXUS"
        throw("$filename is not a Nexus file!")
    end # if

    while true
        line = popfirst!(content)
        if line == "BEGIN DATA;"
            # here comes the important information, so break the loop
            break
        end #if
    end # while
    # get meta info
    ntax, nchar, gap, missing_representation = extract_meta_info(content)
    df = create_nexusdf(content)
    return ntax, nchar, gap, missing_representation, df
end # function ParseNexus

"""
    extract_meta_info(content::Array{String})

This function extracts some meta information from the content of the nexus file.
It "eats-up" the stack
"""
function extract_meta_info(content::Array{String})
    lct::Int64 = 0
    ntax::Int64 = 0
    nchar::Int64 = 0
    gap::String = "-"
    missing_representation::String = "?"
    while true
        line = popfirst!(content)
        if line=="MATRIX"
            # meta info is done, data begins now
            break
        else
            splitted_line = split(line)
            for entry in splitted_line
                info = split(entry, "=")
                if length(info) != 1
                    choped = info[2]
                    if endswith(choped, ";")
                        choped = chop(choped)
                    end #if

                    if info[1] == "ntax"
                        ntax = parse(Int64, choped)
                    elseif info[1] == "NCHAR"
                        nchar = parse(Int64, choped)
                    elseif info[1] == "GAP"
                        gap = choped
                    elseif info[1] == "MISSING"
                        missing_representation = choped
                    else
                        continue
                    end # if
                end # if
            end # for
        end # if
    end # while
    return ntax, nchar, gap, missing_representation
end # function extract_meta_info

"""
    create_nexusdf(filecontent::Array{String})::DataFrame

This function creates a DataFrame of the acutal data.
"""
function create_nexusdf(filecontent::Array{String})::DataFrame
    df = DataFrame(Language=String[], Data=String[])
    while true
        line = popfirst!(filecontent)
        if line == ""
            continue
        elseif line == ";"
            break
        else
            push!(df, split(line))
        end # if
    end # while
    return df
end # function create_nexusdf

"""
    make_tree_with_data(filename::String)::Node

This function creates a tree where the terminal nodes get the data specified in
the NEXUS file.
"""
function make_tree_with_data(filename::String)#::Tree_Module.Node
    # get all the information from the NEXUS file
    n_tax, nc, gap, miss, df = ParseNexus(filename)

    # create random tree
    new_tree = create_tree_from_leaves(df[!,:Language], nc)

    # iterate through the data frame and get the node information
    for row in eachrow(df)
        data_vec = zeros(Float64, (2, nc))
        for (ind, i) in enumerate(row.Data)
            if i == '0'
                data_vec[1,ind] = 1.0
            elseif i == '1'
                data_vec[2,ind] = 1.0
            else
                data_vec[1, ind] = 1.0
                data_vec[2, ind] = 1.0
            end # if
        end # for
        node = find_by_name(new_tree, row.Language)
        node.data = data_vec
    end # for
    return new_tree
end # function make_tree_with_data



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
        if line == "begin data;"
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
        if line=="matrix"
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
                    elseif info[1] == "nchar"
                        nchar = parse(Int64, choped)
                    elseif info[1] == "gap"
                        gap = choped
                    elseif info[1] == "missing"
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

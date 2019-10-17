"""
    ParseNexus(filename::String)

This function parses a NEXUS file which stores the input for the MCMC compuation.
The file should follow the conventions used for MrBayes.
"""
function ParseCSV(filename::String, header::Bool=True)
    open(filename, "r") do file
        global content = readlines(file)
    end # do

    # remove the header if necessary
    header ||  popfirst!(content)

    df = create_csvdf(content)
    dimensions = size(df)
    ntax = dimenseions[1]
    nchar = dimenseions[2]
    return ntax, nchar, gap, missing_representation, df
end # function ParseCSV


"""
    function create_csvdf(filecontnt::Array{String}, separator::AbstractString=",")::DataFrame

This function parses a CSV file and returns its content as a DataFrame.

NOTE: May later on be replaced by the respective DataFrames function.
"""
function create_csvdf(filecontent::Array{String}, separator::AbstractString=",")::DataFrame
    df = DataFrame(Language=String[], Data=String[])
    while !isempty(filecontent)
        line = popfirst!(filecontent)
        push!(df, split(line))
    end # while
    return df
end # function ParseCSV

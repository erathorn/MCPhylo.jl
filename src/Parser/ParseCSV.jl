"""
    ParseCSV(filename::String, header::Bool=true)

This function parses a CSV file containing input for the MCMC compuation.
The file should follow the conventions used for MrBayes. For example:

    Swedish_0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,?,0,0,?,0,0
    Welsh_N_0,0,0,0,0,0,0,?,0,0,0,0,?,?,0,0,?,0,0,0,1,?,?,0
    Sardinian_N_0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,?,0,0,0,0,0

 Set input for "header" to false
if no header is present in the file.
"""
function ParseCSV(filename::String, header::Bool=true)
    open(filename, "r") do file
        global content = readlines(file)
    end # do
    # remove the header if necessary
    header ||  popfirst!(content)
    df = create_csvdf(content)
    dimensions = size(df)
    ntax = dimensions[1]
    nchar = dimensions[2]
    symbols = get_alphabet(df)
    return ntax, nchar, gap, missing_representation, symbols, df
end # function ParseCSV


"""
    function create_csvdf(filecontnt::Array{String}, separator::AbstractString=",")::DataFrame

This function parses a CSV file and returns its content as a DataFrame.

NOTE: May later on be replaced by the respective DataFrames function.
"""
function create_csvdf(filecontent::Array{String}, separator::AbstractString=",")::DataFrame
    language = []
    data = []
    #df's created after editing the filecontent into separate language/data objects
    for line in filecontent
        splitline = split(line,separator)
        push!(language,splitline[1])
        push!(data,splitline[2:end])
    end #for
    df = DataFrame(Language = language,Data = data)
    return df
end # function create_csvdf

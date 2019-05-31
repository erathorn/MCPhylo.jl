module SubstitutionMat

using Markdown

export exponentiate_binary

"""
    exponentiate_binary(pi::Float64, t::float::64, r::Float64)::Array{Float64}

This function returns the expontiatied matrix for the restriction site model.
Following Felsenstein 1981
"""
function exponentiate_binary(pi::Number, t::Float64, r::Float64)::Array{Float64}
    # calculate some auxiliary variables
    ext::Float64 = exp(-t*r)
    ext_::Float64 = 1.0-ext
    p_::Float64 = 1.0-pi
    v_::Float64 = ext_*pi
    w_::Float64 = ext_*p_
    # return the expontiated matrix
    return [ext+v_ w_;
            v_ ext+w_]
end # function

end # moduel SubstitutionMat

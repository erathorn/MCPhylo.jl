module SubstitutionMat

using Markdown

export exponentiate_binary

"""
    exponentiate_binary(pi::Float64, t::float::64, r::Float64)::Array{Float64}

This function returns the expontiatied matrix for the restriction site model.
"""
function exponentiate_binary(pi::Float64, t::Float64, r::Float64)::Array{Float64}
    ext::Float64 = exp(-t*r)
    ext_::Float64 = 1.0-ext
    p_::Float64 = 1.0-pi
    v_::Float64 = ext_*pi
    w_::Float64 = ext_*p_
    return [ext+v_ w_;
            v_ ext+w_]
end # function

end # moduel SubstitutionMat

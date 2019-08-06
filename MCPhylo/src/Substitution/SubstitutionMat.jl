"""
    exponentiate_binary(pi::Float64, t::float::64, r::Float64)::Array{Float64,2}

This function returns the expontiatied matrix for the restriction site model.
Following Felsenstein 1981
"""
function exponentiate_binary(pi_::Number, t::Float64, r::Float64)::Array{Float64,2}
    # calculate some auxiliary variables
    ext::Float64 = exp(-t*r)
    ext_::Float64 = 1.0-ext
    p_::Float64 = 1.0-pi_
    v_::Float64 = ext_*pi_
    w_::Float64 = ext_*p_

    # return the expontiated matrix
    return [ext+v_ w_;
            v_ ext+w_]
end # function

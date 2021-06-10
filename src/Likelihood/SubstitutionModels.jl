const MCP_TIME_MIN = 1.0E-11
const MCP_TIME_MAX = 100.0

"""
    F81(base_freq::Vector{Float64},
                SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}}

Calculate the eigenvalue decomposition of the Q matrix of the F81 model.
The `SubstitutionRates` are ignored, and just for call stability.

The function returns the Eigenvectors, Eigenvalues, inverse of eigenvectors and
    the scale factor for expected number changes per site
"""
function F81(base_freq::Vector{Float64}, SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}, Float64}
    Nbases = length(base_freq)
    Q::Array{Float64,2} = ones(Nbases,Nbases)
    Q[diagind(Nbases,Nbases)] .= 0
    Q .*= reverse(base_freq)
    dia = sum(Q,dims=1)
    Q[diagind(Nbases,Nbases)] = -dia
    Q = transpose(Q)
    D, U = eigen(Q)
    Uinv = inv(U)
    mu::Float64 =  1.0 / (1-dot(base_freq))
    return U, D, Uinv, mu
end

"""
    Restriction(base_freq::Vector{Float64},
                SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}}
Calculate the eigenvalue decomposition of the Q matrix of the restriction site model.
The `SubstitutionRates` are ignored, and just for call stability.

The function returns the Eigenvectors, Eigenvalues, inverse of eigenvectors and
    the scale factor for expected number changes per site
"""
function Restriction(base_freq::Vector{Float64}, SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}, Float64}
    Nbases = length(base_freq)
    Q::Array{Float64,2} = ones(Nbases,Nbases)
    Q[diagind(Nbases,Nbases)] .= -1
    Q .*= transpose(base_freq)

    D, U = eigen(Q)
    Uinv = inv(U)
    mu::Float64 =  1.0 / (2.0 * prod(base_freq))
    return U, D, Uinv, mu
end


"""
    JC(base_freq::Vector{Float64},
        SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}}

Calculate the eigenvalue decomposition of the Q matrix of the Jukes-Cantor model.
The `SubstitutionRates` are ignored, and just for call stability.

The function returns the Eigenvectors, Eigenvalues and inverse of eigenvectors.
"""
function JC(base_freq::Vector{Float64}, SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}, Float64}
    Nbases = length(base_freq)
    Q::Array{Float64,2} = ones(Nbases,Nbases)
    #μ = SubstitutionRates[1]
    off_diag = 1.0/(Nbases)
    diag = off_diag * (Nbases)
    Q .= off_diag
    Q[diagind(Nbases,Nbases)] .= -diag
    D, U = eigen(Q)
    Uinv = inv(U)
    mu = 1/sum(diag)
    return U, D, Uinv, mu
end


"""
    GTR(base_freq::Vector{Float64},
        SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}}

Calculate the eigenvalue decomposition of the Q matrix of the General Time Reversible model.


The function returns the Eigenvectors, Eigenvalues and inverse of eigenvectors.
"""
function GTR(base_freq::Vector{Float64}, SubstitutionRates::Vector{Float64})::Tuple{Array{Float64,2}, Array{Float64,1}, Array{Float64,2}, Float64}
    Nbases = length(base_freq)
    Q::Array{Float64,2} = setmatrix(SubstitutionRates)
    Q = mapslices(x-> x.*base_freq, Q, dims=2)
    dia = sum(Q, dims=2)
    Q[diagind(Nbases,Nbases)] = -dia
    D, U = eigen(Q)
    Uinv = inv(U)
    return U, D, Uinv, 1/sum(dia)
end


"""
    function freeK(base_freq::Vector{Float64},SubstitutionRates::AbstractArray)::Tuple{Array{N,2}, Array{N,1}, Array{N,2}, M} where {N <: Number, M <: Number}

FreeK model of substitution.

"""
function freeK(
      base_freq::Vector{Float64},
      SubstitutionRates::AbstractArray,
      )::Tuple{Array, Array, Array, Float64}# where {N <: Number, M <: Number}
      Nrates = length(SubstitutionRates)
      Nbases = Int(ceil(sqrt(Nrates)))
      Q = zeros(Nbases, Nbases)
      counter = 1
      for i in 1:Nbases, j in 1:Nbases
            if i!=j
                  Q[j,i] = SubstitutionRates[counter]
                  counter += 1
            end
      end
      dia = sum(Q, dims = 2)
      Q[diagind(Nbases, Nbases)] = -dia
      D, U = eigen(Q)
      Uinv = inv(U)
      return U, D, Uinv, 1 / sum(dia)
end


## Helper Functions ##
function setmatrix(vec_vals::Array{T,1})::Array{T,2} where T
    n = length(vec_vals)
    s = round((sqrt(8n+1)+1)/2)
    s*(s-1)/2 == n || error("setmatrix: length of vector is not triangular")
    k = 0
    Q = [ i<j ? (k+=1; vec_vals[k]) : 0 for i=1:s, j=1:s ]
    Q.+transpose(Q)
end


### Calculate Transition Matrices

function calculate_transition(f, rate::R, mu::R, time::R, U::A, Uinv::A, D::Vector, pi_::Vector)::Array{Float64,2} where {R<:Real, A<:AbstractArray{<:Real}}
    return_mat = Array{Float64,2}(undef, length(pi_),length(pi_))
    t = rate * mu * time
    BLAS.gemm!('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* t)), U), Uinv, 0.0, return_mat)
    return return_mat
end

function calculate_transition(f::typeof(JC), rate::R, mu::R, time::R, U::A, Uinv::A, D::Vector, pi_::Vector)::Array{Float64,2} where {R<:Real, A<:AbstractArray{<:Real}}
    return_mat = similar(U)
    t = rate * mu * time
    if t < MCP_TIME_MIN
        return_mat .= 0.0
        return_mat[diagind(return_mat)] .= 1.0
    elseif t > MCP_TIME_MAX
        return_mat .= 1.0/length(pi_)
    else
        BLAS.gemm!('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* t)), U), Uinv, 0.0, return_mat)
    end
    return_mat
end

function calculate_transition(f::typeof(F81), rate::R, mu::R, time::R, U::A, Uinv::A, D::Vector, pi_::Vector)::Array{Float64,2} where {R<:Real, A<:AbstractArray{<:Real}}
    return_mat = similar(U)
    t = rate * mu * time
    if t < MCP_TIME_MIN
        return_mat .= 0.0
        return_mat[diagind(return_mat)] .= 1.0
    elseif t > MCP_TIME_MAX
        return_mat .= reverse(pi_)
        return_mat = collect(transpose(return_mat))
    else
        BLAS.gemm!('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* t)), U), Uinv, 0.0, return_mat)
    end
    return_mat
end



function calculate_transition(f, rate::R, mu::R, time::R, U::A, Uinv::A, D::Vector, pi_::Vector)::Array{Float64,2} where {R<:Real, A<:AbstractArray{<:Complex}}
    return_mat = Array{Float64,2}(undef, length(pi_),length(pi_))
    t = rate * mu * time
    return_mat .= abs.(BLAS.gemm('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* t)), U), Uinv))
    return return_mat
end

function calculate_transition(f::typeof(JC), rate::R, mu::R, time::R, U::A, Uinv::A, D::Vector, pi_::Vector)::Array{Float64,2} where {R<:Real, A<:AbstractArray{<:Complex}}
    return_mat = similar(U)
    t = rate * mu * time
    if t < MCP_TIME_MIN
        return_mat .= 0.0
        return_mat[diagind(return_mat)] .= 1.0
    elseif t > MCP_TIME_MAX
        return_mat .= 1.0/length(pi_)
    else
        return_mat .= abs.(BLAS.gemm('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* t)), U), Uinv))
    end
    return_mat
end

function calculate_transition(f::typeof(Restriction), rate::R, mu::R, time::R, U::A, Uinv::A, D::Vector, pi_::Vector)::Array{Float64,2} where {R<:Real, A<:AbstractArray{<:Complex}}
    return_mat = similar(U)
    t = rate * mu * time
    if t < MCP_TIME_MIN
        return_mat .= 0.0
        return_mat[diagind(return_mat)] .= 1.0
    elseif t > MCP_TIME_MAX
        return_mat .= reverse(pi_)
        return_mat = collect(transpose(return_mat))
    else
        return_mat .= abs.(BLAS.gemm('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* t)), U), Uinv))
    end
    return_mat
end
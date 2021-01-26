

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
    Q .*= reverse(base_freq)
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
    μ = SubstitutionRates[1]
    off_diag = μ/Nbases
    diag = off_diag * (Nbases-1)
    Q .= off_diag
    Q[diagind(Nbases,Nbases)] .= -diag
    println(Q)
    D, U = eigen(Q)
    Uinv = inv(U)
    mu = off_diag[1]
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


function setmatrix(vec_vals::Array{T,1})::Array{T,2} where T
    n = length(vec_vals)
    s = round((sqrt(8n+1)+1)/2)
    s*(s-1)/2 == n || error("setmatrix: length of vector is not triangular")
    k = 0
    Q = [ i<j ? (k+=1; vec_vals[k]) : 0 for i=1:s, j=1:s ]
    Q.+transpose(Q)
end

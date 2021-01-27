
mutable struct PhyloDist <: DiscreteMatrixDistribution
    tree::T where T <: GeneralNode
    base_freq::Array{Float64}
    substitution_rates::Array{Float64}
    rates::Array{Float64}
    substitution_model::Function

    function PhyloDist(tree::T, base_freq::Array{Float64}, substitution_rates::Array{Float64},
                       rates::Array{Float64}, substitution_model::Function) where {T <: GeneralNode}
        new(tree, base_freq, substitution_rates,rates, substitution_model)
    end

end

function PhyloDist(tree::T, base_freq::S, substitution_rates::A, rates::B, substitution_model::Function) where {T<:TreeVariate, S<:DenseArray{Float64}, A<:DenseArray{Float64}, B<:DenseArray{Float64}}
    PhyloDist(tree.value, Array(base_freq), Array(substitution_rates), Array(rates), substitution_model)
end

function PhyloDist(my_tree::T, base_freq::S, substitution_rates::R, rates::R, substitution_model::Function) where {T<:TreeVariate, S<:DenseArray{Float64}, R<:Real}
    PhyloDist(my_tree.value, Array(base_freq), [substitution_rates], [rates], substitution_model)
end





minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

Base.size(d::PhyloDist) = (d.nbase, d.nsites, d.nnodes)

function logpdf(d::PhyloDist, x::AbstractArray)
    mt = post_order(d.tree)
    nba, nsi, nno = size(x)
    data = Array{Float64, 4}(undef, nba, nsi, length(d.rates), nno)
    @inbounds for i in 1:length(d.rates)
        data[:, :, i, :] .= x
    end
    U, D, Uinv, mu = d.substitution_model(d.base_freq, d.substitution_rates)
    FelsensteinFunction(mt, d.base_freq, d.rates, U, D, Uinv, mu, data, false)[1]
end


function gradlogpdf(d::PhyloDist, x::AbstractArray)

    mt = post_order(d.tree)
    nba, nsi, nno = size(x)
    data = Array{Float64, 4}(undef, nba, nsi, length(d.rates), nno)
    @inbounds for i in 1:length(d.rates)
        data[:, :, i, :] .= x
    end
    U, D, Uinv, mu = d.substitution_model(d.base_freq, d.substitution_rates)
    resultate = FelsensteinFunction(mt, d.base_freq, d.rates, U, D, Uinv, mu, data)

    resultate
end

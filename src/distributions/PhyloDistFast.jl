

"""
This structure implements a Distribution whos likelihood is calculated
according to Felsensteins algorithm.
"""
struct PhyloDistFast <: DiscreteMatrixDistribution
    tree::T where {T<:GeneralNode}
    base_freq::Vector{Float64}
    substitution_rates::Vector{Float64}
    rates::Vector{Float64}
    substitution_model::Function
    nbase::Int64
    nnodes::Int64

    function PhyloDistFast(
        tree::T,
        base_freq::A,
        substitution_rates::A,
        rates::A,
        substitution_model::Function,
    ) where {T<:GeneralNode,A<:AbstractArray{<:Real}}
        new(
            tree,
            base_freq,
            substitution_rates,
            rates,
            substitution_model,
            length(base_freq),
            length(post_order(tree)),
        )
    end

end

function PhyloDistFast(
    tree::T,
    base_freq::S,
    substitution_rates::A,
    rates::B,
    substitution_model::Function,
) where {
    T<:TreeVariate,
    S<:DenseArray{Float64},
    A<:DenseArray{Float64},
    B<:DenseArray{Float64},
}
    PhyloDistFast(
        tree.value,
        Array(base_freq),
        Array(substitution_rates),
        Array(rates),
        substitution_model,
    )
end

function PhyloDistFast(
    my_tree::T,
    base_freq::S,
    substitution_rates::R,
    rates::R,
    substitution_model::Function,
) where {T<:TreeVariate,S<:DenseArray{Float64},R<:Real}
    PhyloDistFast(
        my_tree.value,
        Array(base_freq),
        [substitution_rates],
        [rates],
        substitution_model,
    )
end

"""
function PhyloDistFast(my_tree::T, base_freq::S, substitution_rates::R, rates::R, substitution_model::Function) where {T<:GeneralNode, S<:DenseArray{Float64}, R<:Real}

Convenience function which can work with MCPhylo types.
"""
function PhyloDistFast(
    my_tree::T,
    base_freq::S,
    substitution_rates::R,
    rates::R,
    substitution_model::Function,
) where {T<:GeneralNode,S<:DenseArray{Float64},R<:Real}
    PhyloDistFast(
        my_tree,
        Array(base_freq),
        [substitution_rates],
        [rates],
        substitution_model,
    )
end

function PhyloDistFast(
    my_tree::A,
    substitution_rates::B,
    rates::R,
    substitution_model::K,
) where {A<:AbstractNode,B<:DenseArray{Float64},R<:DenseArray{Float64},K<:typeof(freeK)}
    U, D, Uinv, _ = freeK(Float64[], Array(substitution_rates))
    D[:] .= 0
    D[end] = 1
    eq = real.((U*diagm(D)*Uinv)[1, :])
    PhyloDistFast(my_tree, eq, substitution_rates, rates, substitution_model)
end

minimum(d::PhyloDistFast) = -Inf
maximum(d::PhyloDistFast) = Inf

Base.size(d::PhyloDistFast) = (d.nbase, 1, d.nnodes)

function logpdf(d::PhyloDistFast, x::AbstractArray{<:Real,3})::Float64
    mt = post_order(d.tree)
    
    U, D, Uinv, mu = d.substitution_model(
        d.base_freq,
        d.substitution_rates,
    )::Tuple{Matrix,Vector,Matrix,Float64}
    
    #ll = zero(Float64)
    #nsites = size(x, 2)
    #lck = Threads.SpinLock()
    #minparts = min(nsites, 200)
    #parts = max(minparts, Int(round(nsites / Threads.nthreads())))

    #@inbounds for r = 1:length(d.rates)
    #    Threads.@threads for chunk in collect(Iterators.partition(1:nsites, parts))
    ll1 = FelsensteinFunction(
        mt,
        d.base_freq,
        d.rates,
        U,
        D,
        Uinv,
        mu,
        x,
        d.substitution_model,
    )
            
        
    

    ll1
end

function gradlogpdf(d::PhyloDistFast, x::AbstractArray)
    mt = post_order(d.tree)
    U, D, Uinv, mu = d.substitution_model(
        d.base_freq,
        d.substitution_rates,
    )::Tuple{Matrix,Vector,Matrix,Float64}
    ll = zero(Float64)
    gr = zeros(size(x, 3) - 1)
    nsites = size(x, 2)
    lck = Threads.SpinLock()
    minparts = min(nsites, 200)
    parts = max(minparts, Int(round(nsites / Threads.nthreads())))

    @inbounds for r = 1:length(d.rates)
        Threads.@threads for chunk in collect(Iterators.partition(1:nsites, parts))
            ll1, gr1, _ = FelsensteinFunction_wg(
                mt,
                d.base_freq,
                d.rates[r],
                U,
                D,
                Uinv,
                mu,
                x[:, chunk, :],
                d.substitution_model,
            )
            lock(lck) do
                ll += ll1
                gr .+= gr1
            end
        end
    end

    ll, gr
end



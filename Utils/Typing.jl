
Base.convert(::Type{T}, v::TreeVariate) where T<:TreeVariate = v

Base.size(v::Node) = size(post_order(v))



function StochasticTree(f::Function, monitor::Union{Bool, Vector{Int}}=true)
    value = Node()
    fx, src = Mamba.modelfxsrc(Mamba.depfxargs, f)
    s = TreeStochastic(value, :nothing, Int[], fx, src, Symbol[], Mamba.NullUnivariateDistribution())
    #setmonitor!(s, monitor)
end

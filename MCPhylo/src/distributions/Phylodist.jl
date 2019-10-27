
mutable struct PhyloDist <: DiscreteMatrixDistribution
    my_tree
    mypi::Real
    rates::Vector
    nnodes::Int64
    nbase::Int64
    nsites::Int64
end
minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

Base.size(d::PhyloDist) = (d.nnodes, d.nbase, d.nsites)

function logpdf(d::PhyloDist, x::AbstractArray)

    mt = post_order(d.my_tree.value)
    blv = get_branchlength_vector(d.my_tree.value)
    #rates = ones(3132)
    return FelsensteinFunction(mt, d.mypi, d.rates, x, d.nsites, blv)
end

function gradlogpdf(d::PhyloDist, x::AbstractArray)
    #mt = pre_order(d.my_tree.value)
    mt = post_order(d.my_tree.value)
    f(y) = exp(FelsensteinFunction(mt, d.mypi, d.rates, x, d.nsites, y))

    return gradient(f, get_branchlength_vector(d.my_tree.value), :forward)
    #return GradiantLog(mt, d.mypi, d.rates, x, d.nsites)
end

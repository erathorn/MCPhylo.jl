
mutable struct PhyloDist <: DiscreteMatrixDistribution
    my_tree
    mypi::Number
    rates::Vector
    nbase::Int64
    nsites::Int64
    nnodes::Int64
    blv::Vector{Float64}
    function PhyloDist(my_tree, mypi, rates, nbase, nsites, nnodes)
        new(my_tree, mypi, rates, nbase, nsites, nnodes, zeros(nnodes-1))
    end

end
minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

Base.size(d::PhyloDist) = (d.nbase, d.nsites, d.nnodes)

function logpdf(d::PhyloDist, x::AbstractArray)


    mt = post_order(d.my_tree.value)

    get_branchlength_vector(d.my_tree.value, d.blv)
    #get_branchlength_vector(d.my_tree, d.blv)

    return FelsensteinFunction(mt, d.mypi, d.rates, x, d.nsites, d.blv)
end


function gradlogpdf(d::PhyloDist, x::AbstractArray)

    mt = post_order(d.my_tree.value)
    f(y) = FelsensteinFunction(mt, d.mypi, d.rates, x, d.nsites, y)
    get_branchlength_vector(d.my_tree.value, d.blv)
    #gr = gradient(f, d.blv, :forward)

    return val_der(f, d.blv)
end

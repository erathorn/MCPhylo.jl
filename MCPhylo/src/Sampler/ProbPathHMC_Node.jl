function setinits!(d::TreeStochastic, m::Model, x::Node)
    d.value = x
    d.distr = d.eval(m)
    setmonitor!(d, d.monitor)
end # function


#=
tester:
- Julia version: 1.3.0
- Author: erathorn
- Date: 2019-05-07
ohrnstein-uhlenbeck
tree length in compound diriclet
=#

include("./MCPhylo/src/MCPhylo.jl")
using .MCPhylo

using LinearAlgebra
using Random
using StatsFuns
Random.seed!(42)



mt, df = make_tree_with_data("LangData/Dravidian.cc.phy.nex"); # load your own nexus file
#blv = get_branchlength_vector(mt)

#co = MCPhylo.to_covariance(mt, blv)

#cholesky(MCPhylo.to_covariance(mt).*inits[1][:σi][1]^2)
#inits[1][:σi][1]
l = length(MCPhylo.get_leaves(mt))
n = size(df, 2)
arr = [df[:,:,n.num] for n in MCPhylo.get_leaves(mt)]
arrn = zeros(l, n)
for (ind, i) in enumerate(arr)
    arrn[ind, :] .= i[1,:]
end

my_data = Dict{Symbol, Any}(
  :mtree => mt,
  :arrn => arrn,
  :leaves => size(arrn, 1),
  :residuals => size(arrn, 2),
  :nnodes => size(df, 3),
  :σi => ones(3)
);

model =  Model(
    arrn = Stochastic(2,
    (μi, mtree, P, scaler) -> BrownianPhylo(μi, mtree, my_data[:σi], P, scaler,  my_data[:residuals], my_data[:leaves]),false),

    μi = Logical(2, (μ)->(ones(my_data[:residuals],my_data[:leaves]) .= μ)', false),
    μ = Stochastic(1, (μH, σH) -> Normal(μH, σH), false),

    μH = Stochastic(()->Normal()),
    σH = Stochastic(()->Exponential()),
    P = Stochastic(2, ()  -> Normal(),false),

    #σii = Logical(1, σi -> σi.^2,false),
    #σi = Stochastic(1,(λ) -> Exponential(λ), false),
    #σi = Stochastic(1,() -> Gamma(2,2), false),
    scaler = Stochastic(()-> Beta(1,1)),
    #λ = Stochastic(() -> Exponential()),
    mtree = Stochastic(Node(), () -> CompoundDirichlet(1.0,1.0,0.100,1.0), my_data[:nnodes]+1, true),
    #mtree = Stochastic(Node(), () -> MCPhylo.exponentialBL(0.5), my_data[:nnodes]+1, true),

     )


# intial model values
inits = [Dict{Symbol, Union{Any, Real}}(
    :mtree => mt,
    :arrn => my_data[:arrn],
    :nnodes => my_data[:nnodes],
    :nl => my_data[:nnodes],
    :nsites => my_data[:residuals],
    :P => randn(my_data[:leaves], my_data[:residuals]),
    :μH => rand(),
    :σH => rand(),
    :μ => randn(my_data[:residuals]),
    :λ => 1,
    :σi => ones(my_data[:residuals]),
    :scaler => 0.02#one(Float64)#rand()
    )
    ]

scheme = [#RWM(:mtree, 1),
          #Slice(:mtree, 0.05, Multivariate),
          PNUTS(:mtree),
          #Slice(:μH, 0.05, Univariate),

          Slice(:P, 0.1, Multivariate),
          #RWM(:P, 1),
          #SliceSimplex(:ζ),
          #, :λ
          #Slice(:σi, 0.05, Univariate),
          #:σH,
          Slice([:σH,:μH ,:scaler], 0.05, Multivariate)
          ]

setsamplers!(model, scheme);
using Zygote

Zygote.@nograd isposdef
# do the mcmc simmulation. if trees=true the trees are stored and can later be
# flushed ot a file output.
sim = mcmc(model, my_data, inits, 20, burnin=10,thin=2, chains=1, trees=true)
#sim = mcmc(sim, 10, trees=true)

# write the output to a path specified as the second argument
to_file(sim, "multiv", 2)


m = (ones(my_data[:leaves],my_data[:residuals]) .= transpose(inits[1][:μ]))
blv = get_branchlength_vector(mt)
cov_fun_mat = MCPhylo.to_covariance_func(mt)
si = inits[1][:σi]
P = inits[1][:P]
f(y) = mlpdf_h(Array(m), cov_fun_mat, y, si, P, 0.02, my_data[:residuals], arrn)

r = Zygote.pullback(tester1, cov_fun_mat, blv, Array(m), P, 231, 38,0.02, arrn)
r[2](1.0)
using BenchmarkTools

function test_grad(cov_fun_mat, blv, m, P, nres, nchar,scaler, arrn)
    f(y) = tester1(cov_fun_mat, y, Array(m), P, nres, nchar,scaler, arrn)
    r = Zygote.pullback(f, blv)
    r[1], r[2](1.0)[1]
end
r = Zygote.pullback(tester1, cov_fun_mat, blv, Array(m), P, 231, 38,0.02, arrn)
@profiler r[2](1.0)
@benchmark tester1(cov_fun_mat, blv, m, P, 231, 38,0.02, arrn)
@benchmark test_grad(cov_fun_mat, blv, m, P, 231, 38,0.02, arrn)
test_grad(cov_fun_mat, blv, m, P, 231, 38,0.02, arrn)
function p_apply_array(f::Function, x::Array{T})::T where T<:Real
 f(x)
end
function tester1(cov_fun_mat::Array{Function,2}, blv::Vector{Float64}, mu::Array{Float64,2},
            P::Array{Float64,2}, n_res::Int64, n_langs::Int64, scaler::Float64, data::Array{Float64,2})

    mycov = p_apply_array.(cov_fun_mat, Ref(blv))
    ch = LinearAlgebra.cholesky(mycov)

    chl = ch.L

    ichl = inv(chl)
    si = transpose(ichl)*ichl
    ds = logdet(chl)
    v1 = -n_langs*MCPhylo.l2pi-ds
    result = zero(Float64)
    r = similar(blv)
    @views @inbounds for i=1:n_res
        r = chl * P[:,i]
        result = result + sum(MCPhylo.mymvlogpdf_2(v1, si, r))::Float64
        for j in 1:n_langs
            result = result + my_bernoulli_logpdf_h(myloginvlogit(r[j] + mu[j], scaler), data[j,i])::Float64
        end
    end
    result
end

@inline myinvlogit(x::Real, λ::Real=1.0) = 1.0 / (exp(-λ*x) + 1.0)
@inline function my_bernoulli_logpdf_h(θ::T, x::S)::T where {S <: Real, T<: Real}
    x == 1 ? θ : log(1-exp(θ))
end


@inline myloginvlogit_2(x::Real, λ::Real=1.0) = -log(exp(-λ*x)+1.0)

Zygote.@adjoint function myinvlogit(x::Real, λ::Real)
    myinvlogit(x, λ), Δ -> (nothing, invlogder(Δ, λ))
end

Zygote.@adjoint function myloginvlogit(x::Real, λ::Real)
    myloginvlogit(x, λ), Δ -> (nothing, loginvlogitder(Δ, λ))
end


function loginvlogitder(x::Real, λ::Real)
    λ/(exp(λ*x)+1)
end

function invlogder(x::T, λ::T) where T <: Real
    invlogit(x, λ)*invlogit(-x, λ)
end

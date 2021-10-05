#### Structs for positions in the hamiltonian space

abstract type HMC_State end

mutable struct Tree_HMC_State{T<:GeneralNode} <: HMC_State
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
end


mutable struct Array_HMC_State{T<:Array{<:Real}} <: HMC_State
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
end

function transfer(s1::T)::T where T<:HMC_State
    T(deepcopy(s1.x), s1.r[:],s1.g[:],s1.lf)
end

function transfer(s1::T)::T where T<:Array_HMC_State
    T(s1.x[:], s1.r[:],s1.g[:],s1.lf)
end

hamiltonian(s::HMC_State) = s.lf - 0.5 * dot(s.r)

function transfer!(s1::HMC_State, s2::HMC_State)
    s1.x = deepcopy(s2.x)
    s1.r[:] .= s2.r[:]
    s1.g[:] .= s2.g[:]
    s1.lf = s2.lf
    nothing
end

function transfer!(s1::Array_HMC_State, s2::Array_HMC_State)
    s1.x[:] .= s2.x[:]
    s1.r[:] .= s2.r[:]
    s1.g[:] .= s2.g[:]
    s1.lf = s2.lf
    nothing
end


#### PNUTS Parameter Structs

struct NUTS_StepParams
    μ::Float64
    δ::Float64
    γ::Float64
    κ::Float64
    t0::Int
    τ::Float64
end

function update_step(p::NUTS_StepParams, μ::Float64)
    p1 = NUTS_StepParams(μ, p.δ, p.γ, p.κ, p.t0, p.τ)
    p1
end

mutable struct NUTSstepadapter
    m::Int
    s_bar::Float64
    x_bar::Float64
    params::NUTS_StepParams
    metro_acc_prob::Float64
    avg_nni::Float64

    function NUTSstepadapter(m, s_bar, x_bar, μ, δ, γ, κ, t0, τ)
        new(m, s_bar, x_bar, NUTS_StepParams(μ, δ, γ, κ, t0, τ), 0.0, 0.0)
    end

    function NUTSstepadapter(m, s_bar, x_bar, params::NUTS_StepParams)
        new(m, s_bar, x_bar, params, 0.0, 0.0)
    end
end


mutable struct NUTSMeta
    nni::Float64
    alpha::Float64
    nalpha::Float64
    accnni::Float64
end

NUTSMeta() = NUTSMeta(0.0,0.0,0.0,0.0)

function update!(x::NUTSMeta, y::NUTSMeta)
    x.nni += y.nni
    x.alpha += y.alpha
    x.nalpha += y.nalpha
    x.accnni += y.accnni
end

#### Structs for positions in the hamiltonian space

mutable struct Tree_HMC_State{T<:GeneralNode}
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
end

function transfer(s1::Tree_HMC_State)
    Tree_HMC_State(deepcopy(s1.x), s1.r[:],s1.g[:],s1.lf)
end


hamiltonian(s::Tree_HMC_State) = s.lf - 0.5 * dot(s.r)

function transfer!(s1::Tree_HMC_State, s2::Tree_HMC_State)
    s1.x = deepcopy(s2.x)
    s1.r[:] .= s2.r[:]
    s1.g[:] .= s2.g[:]
    s1.lf = s2.lf
    nothing
end


#### PNUTS Parameter Structs

struct PNUTS_StepParams
    μ::Float64
    δ::Float64
    γ::Float64
    κ::Float64
    t0::Int
    τ::Float64
end

function update_step(p::PNUTS_StepParams, μ::Float64)
    p1 = PNUTS_StepParams(μ, p.δ, p.γ, p.κ, p.t0, p.τ)
    p1
end

mutable struct PNUTSstepadapter
    m::Int
    s_bar::Float64
    x_bar::Float64
    params::PNUTS_StepParams
    metro_acc_prob::Float64
    avg_nni::Float64

    function PNUTSstepadapter(m, s_bar, x_bar, μ, δ, γ, κ, t0, τ)
        new(m, s_bar, x_bar, PNUTS_StepParams(μ, δ, γ, κ, t0, τ), 0.0, 0.0)
    end

    function PNUTSstepadapter(m, s_bar, x_bar, params::PNUTS_StepParams)
        new(m, s_bar, x_bar, params, 0.0, 0.0)
    end
end


mutable struct PNUTSMeta
    nni::Float64
    alpha::Float64
    nalpha::Float64
    accnni::Float64
end

PNUTSMeta() = PNUTSMeta(0.0,0.0,0.0,0.0)

function update!(x::PNUTSMeta, y::PNUTSMeta)
    x.nni += y.nni
    x.alpha += y.alpha
    x.nalpha += y.nalpha
    x.accnni += y.accnni
end

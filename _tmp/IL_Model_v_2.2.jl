begin
    using ModelingToolkit, MethodOfLines, DomainSets, NonlinearSolve
    using GarishPrint
    using Plots
    using LinearAlgebra, SparseArrays, LinearSolve, SparseDiffTools
    using StaticArrays
    using DataInterpolations
    using CoolProp
    using BenchmarkTools
    using ProgressMeter
    using ForwardDiff
    using Optim
    # include("props.jl")
    include("props_regression.jl")
    using Parameters: @unpack
    using ComponentArrays
end
# =======================================================================================
begin "constans"
    #= 
    Nellis, G., & Klein, S. (2008). Heat Transfer. Cambridge: Cambridge University Press.
    EXAMPLE 9.2-1: Diffusion Coefficient for Air-Water Vapor Mixtures
    Bolz, R.E. and G.L. Tuve, Handbook of Tables for Applied Engineering Science, 2nd edition,CRC Press, (1976).
    =#
    Dₐ(T) = -2.775e-6 + 4.479e-8 * T + 1.656e-10 * T^2
    #= 
    Nellis, G., & Klein, S. (2008). Heat Transfer. Cambridge: Cambridge University Press.
    Infinite Dilution Diffusion Coefficients for Liquids
    A modified form of the Tyn-Calus correlation  Eq. 9-31
    Poling, B.E., J.M. Prausnitz, and J. O' Connell, The Properties of Gases and Liquids, 5th Edition,
    McGraw-Hill, New York, (2000), ISBN 0070116822 / 9780070116825. Eq. (11-9.4)
    =#
    ρ_water(T) = CoolProp.PropsSI("D", "T", T, "P", 101325.0, "Water") 
    # https://www.engineeringtoolbox.com/water-surface-tension-d_597.html
    σ_water(T) = (-1e-05 * T^2 - 0.0121 * T + 11.655) * 0.01
    v_water(T) = 18.01528 / 1000ρ_water(T)
    const IL_MW_Base = 25.0
    IL_MW(ξ) = 18.01528 * (1 - ξ) + IL_MW_Base * ξ
    vₗ(T ,ξ) = IL_MW_Base / 1000.0_ρₛₒₗ(T ,ξ)
    # Dₗ(T ,ξ) = 9.013e-16 * (v_water(T) ^ 0.267 / vₗ(T ,ξ) ^ 0.433) * (T / μₛₒₗ(T ,ξ)) * (σₛₒₗ(T ,ξ) / σ_water(T)) ^ 0.15
    Dₗ(T,ξ) = 9.013e-16 * (v_water(T) ^ 0.267 / vₗ(T ,ξ) ^ 0.433) * (T / _μₛₒₗ(T ,ξ)) * (_σₛₒₗ(T ,ξ) / σ_water(T)) ^ 0.15

    const T_w = 15.86 + 273.15 # Evaporator wall temperature
    const ΔT_supersub = 0.0 # Subcooling temperature
    # const T_w = 40.0 + 273.15 # condenser wall temperature
    const N_fin = 48
    # const N_fin = 75
    const MR = 0.052959106 / 0.025182933
    # const MR = 0.063437 / 0.037148
    # const ṁₐ_ₜₒₜ = 0.8 * 0.018709069
    const ṁₐ_ₜₒₜ = 0.025182933
    # const ṁₐ_ₜₒₜ = 0.037148
    const ṁₐ = (ṁₐ_ₜₒₜ / N_fin) * 0.5 # mass flow rate for half of the fin space
    const ṁₛₒₗ = ṁₐ * MR
    const FD = 0.205
    const Tₛₒₗ_ᵢₙ = 22.38 + 273.15
    const ξₛₒₗ_ᵢₙ = 0.28
    const ρₛₒₗ = _ρₛₒₗ(Tₛₒₗ_ᵢₙ, 1 - ξₛₒₗ_ᵢₙ)
    const g = 9.81
    const μₛₒₗ = _μₛₒₗ(T_w, 1 - ξₛₒₗ_ᵢₙ)
    const νₛₒₗ = μₛₒₗ / ρₛₒₗ
    const δₛₒₗ = ∛(3 * ṁₛₒₗ * νₛₒₗ / (ρₛₒₗ * g * FD))
    const H = 0.132
    const FS = 0.00254
    const Uₛₒₗ_ᵣ = ṁₛₒₗ / (ρₛₒₗ * δₛₒₗ * FD)
    const ARₛₒₗ = δₛₒₗ / H
    const Reₛₒₗ = Uₛₒₗ_ᵣ * δₛₒₗ / νₛₒₗ
    const 𝑘ₛₒₗ = _𝑘ₛₒₗ(T_w, 1 - ξₛₒₗ_ᵢₙ)
    const cpₛₒₗ = _cpₛₒₗ(T_w, 1 - ξₛₒₗ_ᵢₙ)
    const Prₛₒₗ = 1.1 * cpₛₒₗ * μₛₒₗ / 𝑘ₛₒₗ 
    const Dₛₒₗ = Dₗ(0.5 * (Tₛₒₗ_ᵢₙ + T_w) ,1 - ξₛₒₗ_ᵢₙ)
    const Scₛₒₗ = 2_000.0
    # const Scₛₒₗ = νₛₒₗ / Dₛₒₗ

    const Tₐ_ᵢₙ = 28.0 + 273.15
    const ωₐ_ᵢₙ = 0.019491
    const ρₐ = _ρₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    const μₐ = _μₐ(Tₐ_ᵢₙ)
    const νₐ = μₐ / ρₐ
    const 𝑘ₐ = _kₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    const αₐ = _αₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    const Prₐ = νₐ / αₐ
    const Scₐ = νₐ / Dₐ(Tₐ_ᵢₙ)
    const δₐ = 0.5 * FS - δₛₒₗ
    const Uₐ_ᵣ = ṁₐ / (ρₐ * δₐ * FD)
    const Reₐ = Uₐ_ᵣ * δₐ / νₐ
    const ARₐ = δₐ / H
    const uᵢₙₜ = 0.5g * δₛₒₗ^2 / νₛₒₗ
    const dpdx = -(3.0 * μₐ * uᵢₙₜ / (δₐ^2)) - (3.0 * μₐ * ṁₐ / (ρₐ * (δₐ ^ 3) * FD))
    const ΔTₐ_ᵣ = Tₐ_ᵢₙ - T_w
    const ΔTₛₒₗ_ᵣ = Tₛₒₗ_ᵢₙ - T_w
    const coeff_ωꜛₐ_ᵢₙₜ = Scₐ * ARₐ * μₛₒₗ * ξₛₒₗ_ᵢₙ / (Scₛₒₗ * ARₛₒₗ * μₐ * ωₐ_ᵢₙ)
    const coeff_Θꜛₐ_ᵢₙₜ_₁ = 𝑘ₛₒₗ * ΔTₛₒₗ_ᵣ * ARₐ / (𝑘ₐ * ΔTₐ_ᵣ * ARₛₒₗ)
    _coeff_Θꜛₐ_ᵢₙₜ_2(Tᵢₙₜ) = (Prₛₒₗ / Scₛₒₗ) * i_fg(Tᵢₙₜ) * ξₛₒₗ_ᵢₙ / (ΔTₛₒₗ_ᵣ * cpₛₒₗ)
end
# 𝑘ₛₒₗ_ = _𝑘ₛₒₗ(T_w, 1 - ξₛₒₗ_ᵢₙ)
# cpₛₒₗ_ = _cpₛₒₗ(T_w, 1 - ξₛₒₗ_ᵢₙ)
# μₛₒₗ_ = _μₛₒₗ(T_w, 1 - ξₛₒₗ_ᵢₙ)
# Prₛₒₗ_ = cpₛₒₗ_ * μₛₒₗ_ / 𝑘ₛₒₗ_ 

# 𝑘ₛₒₗ_ = _𝑘ₛₒₗ(Tₛₒₗ_ᵢₙ, 1 - ξₛₒₗ_ᵢₙ)
# cpₛₒₗ_ = _cpₛₒₗ(Tₛₒₗ_ᵢₙ, 1 - ξₛₒₗ_ᵢₙ)
# μₛₒₗ_ = _μₛₒₗ(Tₛₒₗ_ᵢₙ, 1 - ξₛₒₗ_ᵢₙ)
# Prₛₒₗ_ = cpₛₒₗ_ * μₛₒₗ_ / 𝑘ₛₒₗ_ 


Prₛₒₗ
δₛₒₗ
Reₛₒₗ
Reₐ
ṁₐ
ρₛₒₗ
coeff_ωꜛₐ_ᵢₙₜ
coeff_Θꜛₐ_ᵢₙₜ_₁
_coeff_Θꜛₐ_ᵢₙₜ_2(290)
 0.062513634 / 0.018709069
 νₛₒₗ / Dₛₒₗ
@show Pe1 = Reₛₒₗ * Prₛₒₗ
@show Pe_mass = Reₛₒₗ * Scₛₒₗ
@show Pe2 = Reₐ * Prₐ
@show Pe_mass = Reₐ * Scₐ
# =======================================================================================
begin "kernel"
    function Uꜛₛₒₗ(Yꜛₛₒₗ)
        y = Yꜛₛₒₗ * δₛₒₗ
        uₛₒₗ = g * y * (δₛₒₗ - 0.5 * y) / νₛₒₗ
        return uₛₒₗ / Uₛₒₗ_ᵣ
    end

    function Uꜛₐ(Yꜛₐ)
        y = Yꜛₐ * δₐ
        uₐ = -uᵢₙₜ - (0.5 / μₐ) * dpdx *(δₐ ^ 2 - y ^ 2)
        return uₐ / Uₐ_ᵣ
    end

    const γ = 1.0
    const Δx = 0.01
    const Δy = 0.01
    const M = Int(1.0 / Δx)  # 20 nodes - Forward nodes 
    const N = Int(1.0 / Δy) - 1  # 19 nodes - Intermediate nodes
    const K1 = γ * Δx / Δy^2

    const β_Θₐ = K1 / (Reₐ * Prₐ * ARₐ)
    const β_ωₐ = K1 / (Reₐ * Scₐ * ARₐ)
    const β_Θₛₒₗ = K1 / (Reₛₒₗ * Prₛₒₗ * ARₛₒₗ)
    const β_ξₛₒₗ = K1 / (Reₛₒₗ * Scₛₒₗ * ARₛₒₗ)

    const Θꜛₛₒₗ_0 = collect(0:Δx:1.0) .* (ΔT_supersub / ΔTₛₒₗ_ᵣ)

    function TDMA!(Dl,D,Du,B,X,N)
        @inbounds @simd for i in 2:N
            ω = Dl[i] / D[i-1]
            D[i] = D[i] - ω * Du[i-1]
            B[i] = B[i] - ω * B[i-1]
        end
        @inbounds X[N] = B[N] / D[N]
        @inbounds @simd for i in N-1:-1:1
            X[i] = (B[i] - Du[i] * X[i+1]) / D[i]
        end
    end

    abstract type BCType end
    struct Dirichlet <: BCType end
    struct Neumann <: BCType end

    abstract type BCLocation end
    struct StartLoc <: BCLocation end
    struct EndLoc <: BCLocation end

    # --------------------------------------------------------------------------------------------
    function BC!(Dl,D,Du,B,bc_val,Δ,::Neumann,::StartLoc)
        @inbounds D[1] = D[1] + 4.0 * Dl[1] / 3.0
        @inbounds Du[1] = Du[1] - Dl[1] / 3.0
        @inbounds B[1] = B[1] + Dl[1] * bc_val * Δ * 2.0 / 3.0
        nothing
    end

    function BC!(Dl,D,Du,B,bc_val,Δ,N::Int,::Neumann,::EndLoc)
        @inbounds B[N] = B[N] - Du[N] * bc_val * Δ * 2.0 / 3.0
        @inbounds D[N] = D[N] + 4.0 * Du[N] / 3.0
        @inbounds Dl[N] = Dl[N] - Du[N] / 3.0
        nothing
    end

    function BC!(Dl,D,Du,B,bc_val,Δ,::Dirichlet,::StartLoc) 
        @inbounds B[1] = B[1] - Dl[1] * bc_val
        nothing
    end

    function BC!(Dl,D,Du,B,bc_val,Δ,N::Int,::Dirichlet,::EndLoc) 
        @inbounds B[N] = B[N] - Du[N] * bc_val
        nothing
    end
    # --------------------------------------------------------------------------------------------

    function BC_Numan!(X,bc_val,Δ,::Neumann,::StartLoc) 
        @inbounds X[1] = -2.0 * Δ * bc_val / 3.0 + 4 * X[2] / 3.0 - X[3] / 3.0
        nothing
    end

    function BC_Numan!(X,bc_val,Δ,N::Int,::Neumann,::EndLoc) 
        @inbounds X[N+2] = 2 * Δ * bc_val / 3.0 + 4 * X[N+1] / 3.0 - X[N] / 3.0
        nothing
    end

    function BC_Dirichlet!(X,bc_val,Δ,N::Int,::Dirichlet,::EndLoc) 
        @inbounds X[N+2] = bc_val
        nothing
    end

    function BC_Dirichlet!(X,bc_val,Δ,::Dirichlet,::StartLoc) 
        @inbounds X[1] = bc_val
        nothing
    end

end
# --------------------------------------------------------------------------------------------

function _solve_ξ!(ξꜛₛₒₗ,ξꜛₛₒₗ_ᵢₙₜ,B,α_yₛₒₗ)
    # ξꜛₛₒₗ_ᵢₙₜ = view(X,1:M)
    ∂ξꜛₛₒₗ_∂Yꜛₛₒₗ_0 = 0.0 # ∂ξₛₒₗ / ∂Yₛₒₗ at y = 0
    for i in 1:M
        Du_ξₛₒₗ = [β_ξₛₒₗ / α_yₛₒₗ[k] for k in 1:N] # Upper diagonal of ξₛₒₗ
        Dl_ξₛₒₗ = copy(Du_ξₛₒₗ) # Lower diagonal of ξₛₒₗ
        D_ξₛₒₗ = [-1.0 - 2β_ξₛₒₗ / α_yₛₒₗ[p] for p in 1:N] # Diagonal of ξₛₒₗ
        B_ξₛₒₗ = copy(B) # Right hand side of ξₛₒₗ
        @inbounds @simd for j in 1:N
            B_ξₛₒₗ[j] = -view(ξꜛₛₒₗ,1:N+2,i)[j+1] - β_ξₛₒₗ * ((1.0 - γ) / γ) * (view(ξꜛₛₒₗ,1:N+2,i)[j+2] - 2 * view(ξꜛₛₒₗ,1:N+2,i)[j+1] + view(ξꜛₛₒₗ,1:N+2,i)[j]) / α_yₛₒₗ[j]
        end
        @inbounds BC!(Dl_ξₛₒₗ,D_ξₛₒₗ,Du_ξₛₒₗ,B_ξₛₒₗ,ξꜛₛₒₗ_ᵢₙₜ[i],Δy,N,Dirichlet(),EndLoc()) # Neumann BC at y = 1.0
        @inbounds BC!(Dl_ξₛₒₗ,D_ξₛₒₗ,Du_ξₛₒₗ,B_ξₛₒₗ,∂ξꜛₛₒₗ_∂Yꜛₛₒₗ_0,Δy,Neumann(),StartLoc()) # Neumann BC at y = 0.0
        @inbounds BC_Dirichlet!(view(ξꜛₛₒₗ,:,i+1),ξꜛₛₒₗ_ᵢₙₜ[i],Δy,N,Dirichlet(),EndLoc())
        TDMA!(Dl_ξₛₒₗ,D_ξₛₒₗ,Du_ξₛₒₗ,B_ξₛₒₗ,view(ξꜛₛₒₗ,2:N+1,i+1),N)
        @inbounds BC_Numan!(view(ξꜛₛₒₗ,:,i+1),∂ξꜛₛₒₗ_∂Yꜛₛₒₗ_0,Δy,Neumann(),StartLoc())
    end
end

function _solve_Θ_sol!(Θꜛₛₒₗ,Θꜛₛₒₗ_ᵢₙₜ,α_yₛₒₗ)

    for i in 1:M
        Du_Θₛₒₗ = [β_Θₛₒₗ / α_yₛₒₗ[k] for k in 1:N] # Upper diagonal of Θₛₒₗ
        Dl_Θₛₒₗ = copy(Du_Θₛₒₗ) # Lower diagonal of Θₛₒₗ
        D_Θₛₒₗ = [-1.0 - 2β_Θₛₒₗ / α_yₛₒₗ[p] for p in 1:N] # Diagonal of Θₛₒₗ
        B_Θₛₒₗ = copy(B) # Right hand side of Θₛₒₗ
        @inbounds @simd for j in 1:N
            B_Θₛₒₗ[j] = -view(Θꜛₛₒₗ,1:N+2,i)[j+1] - β_Θₛₒₗ * ((1.0 - γ) / γ) * (view(Θꜛₛₒₗ,1:N+2,i)[j+2] - 2 * view(Θꜛₛₒₗ,1:N+2,i)[j+1] + view(Θꜛₛₒₗ,1:N+2,i)[j]) / α_yₛₒₗ[j]
        end
        @inbounds BC!(Dl_Θₛₒₗ,D_Θₛₒₗ,Du_Θₛₒₗ,B_Θₛₒₗ,Θꜛₛₒₗ_ᵢₙₜ[i],Δy,N,Dirichlet(),EndLoc()) # Dirichlet BC at y = 1.0
        @inbounds BC!(Dl_Θₛₒₗ,D_Θₛₒₗ,Du_Θₛₒₗ,B_Θₛₒₗ,Θꜛₛₒₗ_0[i+1],Δy,Dirichlet(),StartLoc()) # Dirichlet BC at y = 0.0
        @inbounds BC_Dirichlet!(view(Θꜛₛₒₗ,:,i+1),Θꜛₛₒₗ_ᵢₙₜ[i],Δy,N,Dirichlet(),EndLoc()) # Dirichlet BC at y = 1.0
        @inbounds BC_Dirichlet!(view(Θꜛₛₒₗ,:,i+1),Θꜛₛₒₗ_0[i+1],Δy,Dirichlet(),StartLoc()) # Dirichlet BC at y = 0.0
        TDMA!(Dl_Θₛₒₗ,D_Θₛₒₗ,Du_Θₛₒₗ,B_Θₛₒₗ,view(Θꜛₛₒₗ,2:N+1,i+1),N)
    end
end

function _solve_Θ_a!(Θꜛₐ,Θꜛₐ_ᵢₙₜ,α_yₐ)
    ∂Θꜛₐ_∂Yꜛₐ = 0.0 # ∂Θₐ / ∂Yₐ at interface nodes y = 0.0 or Yꜛₐ = 0.0
    for i in 1:M
        Du_Θₐ = [β_Θₐ / α_yₐ[k] for k in 1:N] # Upper diagonal of Θₐ
        Dl_Θₐ = copy(Du_Θₐ) # Lower diagonal of Θₐ
        D_Θₐ = [-1.0 - 2β_Θₐ / α_yₐ[p] for p in 1:N] # Diagonal of Θₐ
        B_Θₐ = copy(B) # Right hand side of Θₐ
        @inbounds @simd for j in 1:N
            B_Θₐ[j] = -view(Θꜛₐ,1:N+2,i)[j+1] - β_Θₐ * ((1.0 - γ) / γ) * (view(Θꜛₐ,1:N+2,i)[j+2] - 2 * view(Θꜛₐ,1:N+2,i)[j+1] + view(Θꜛₐ,1:N+2,i)[j]) / α_yₐ[j]
        end
        @inbounds BC!(Dl_Θₐ,D_Θₐ,Du_Θₐ,B_Θₐ,Θꜛₐ_ᵢₙₜ[i],Δy,N,Dirichlet(),EndLoc()) # Dirichlet BC at y = 1.0
        @inbounds BC!(Dl_Θₐ,D_Θₐ,Du_Θₐ,B_Θₐ,∂Θꜛₐ_∂Yꜛₐ,Δy,Neumann(),StartLoc()) # Neumann BC at y = 0.0
        @inbounds BC_Dirichlet!(view(Θꜛₐ,:,i+1),Θꜛₐ_ᵢₙₜ[i],Δy,N,Dirichlet(),EndLoc()) # Dirichlet BC at y = 1.0
        TDMA!(Dl_Θₐ,D_Θₐ,Du_Θₐ,B_Θₐ,view(Θꜛₐ,2:N+1,i+1),N)
        @inbounds BC_Numan!(view(Θꜛₐ,:,i+1),∂Θꜛₐ_∂Yꜛₐ,Δy,Neumann(),StartLoc()) # Neumann BC at y = 0.0
    end
end

function _solve_ω!(ωꜛₐ,ωꜛ_int,α_yₐ)
    ∂ωꜛₐ_∂Yꜛₐ = 0.0 # ∂ωₐ / ∂Yₐ at interface nodes y = 0.0 or Yꜛₐ = 0.0
    for i in 1:M
        Du_ωₐ = [β_ωₐ / α_yₐ[k] for k in 1:N] # Upper diagonal of ωₐ
        Dl_ωₐ = copy(Du_ωₐ) # Lower diagonal of ωₐ
        D_ωₐ = [-1.0 - 2β_ωₐ / α_yₐ[p] for p in 1:N] # Diagonal of ωₐ
        B_ωₐ = copy(B) # Right hand side of ωₐ
        @inbounds @simd for j in 1:N
            B_ωₐ[j] = -view(ωꜛₐ,1:N+2,i)[j+1] - β_ωₐ * ((1.0 - γ) / γ) * (view(ωꜛₐ,1:N+2,i)[j+2] - 2 * view(ωꜛₐ,1:N+2,i)[j+1] + view(ωꜛₐ,1:N+2,i)[j]) / α_yₐ[j]
        end
        @inbounds BC!(Dl_ωₐ,D_ωₐ,Du_ωₐ,B_ωₐ,ωꜛ_int[i],Δy,N,Dirichlet(),EndLoc()) # Dirichlet BC at y = 1.0
        @inbounds BC!(Dl_ωₐ,D_ωₐ,Du_ωₐ,B_ωₐ,∂ωꜛₐ_∂Yꜛₐ,Δy,Neumann(),StartLoc()) # Neumann BC at y = 0.0
        @inbounds BC_Dirichlet!(view(ωꜛₐ,:,i+1),ωꜛ_int[i],Δy,N,Dirichlet(),EndLoc()) # Dirichlet BC at y = 1.0
        TDMA!(Dl_ωₐ,D_ωₐ,Du_ωₐ,B_ωₐ,view(ωꜛₐ,2:N+1,i+1),N)
        @inbounds BC_Numan!(view(ωꜛₐ,:,i+1),∂ωꜛₐ_∂Yꜛₐ,Δy,Neumann(),StartLoc()) # Neumann BC at y = 0.0
    end
end

function solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int,B,α_yₛₒₗ,α_yₐ)
    ξꜛₛₒₗ_ᵢₙₜ = view(X,1:M)
    Θꜛₛₒₗ_ᵢₙₜ = view(X,M+1:2M)
    _solve_ξ!(ξꜛₛₒₗ,ξꜛₛₒₗ_ᵢₙₜ,B,α_yₛₒₗ)
    _solve_Θ_sol!(Θꜛₛₒₗ,Θꜛₛₒₗ_ᵢₙₜ,α_yₛₒₗ)
    # --------------------------------------------------------------------------------------------
    ξ_int .= ξꜛₛₒₗ_ᵢₙₜ .* ξₛₒₗ_ᵢₙ
    T_int .= Θꜛₛₒₗ_ᵢₙₜ .* ΔTₛₒₗ_ᵣ .+ T_w
    Θꜛₐ_ᵢₙₜ .= (T_int .- T_w) ./ ΔTₐ_ᵣ |> reverse!
    ωꜛ_int .= (0.62185 .* _Pᵥₐₚₒᵣ_ₛₒₗ.(T_int,1.0 .- ξ_int) ./ (101325.0 .- _Pᵥₐₚₒᵣ_ₛₒₗ.(T_int,1.0 .- ξ_int))) ./ ωₐ_ᵢₙ |> reverse!
    # --------------------------------------------------------------------------------------------
    _solve_Θ_a!(Θꜛₐ,Θꜛₐ_ᵢₙₜ,α_yₐ)
    _solve_ω!(ωꜛₐ,ωꜛ_int,α_yₐ)    
    # --------------------------------------------------------------------------------------------
    first_der_at_int = (x , x⁻ , x⁻⁻) -> (3x - 4x⁻ + x⁻⁻) / 2Δy
    ∂ξₛₒₗ_∂Yₛₒₗ_int = map(first_der_at_int,ξꜛₛₒₗ[N+2,2:M+1],ξꜛₛₒₗ[N+1,2:M+1],ξꜛₛₒₗ[N,2:M+1])
    ∂Θₛₒₗ_∂Yₛₒₗ_int = map(first_der_at_int,Θꜛₛₒₗ[N+2,2:M+1],Θꜛₛₒₗ[N+1,2:M+1],Θꜛₛₒₗ[N,2:M+1])
    ∂Θₐ_∂Yₐ_int = map(first_der_at_int,Θꜛₐ[N+2,2:M+1],Θꜛₐ[N+1,2:M+1],Θꜛₐ[N,2:M+1]) |> reverse
    ∂ωₐ_∂Yₐ_int = map(first_der_at_int,ωꜛₐ[N+2,2:M+1],ωꜛₐ[N+1,2:M+1],ωꜛₐ[N,2:M+1]) |> reverse
    loss1 = map((x,y) -> x + coeff_ωꜛₐ_ᵢₙₜ * y, ∂ωₐ_∂Yₐ_int[1:M],∂ξₛₒₗ_∂Yₛₒₗ_int[1:M])
    loss2 = map((x,y,T,z) -> begin
                        - coeff_Θꜛₐ_ᵢₙₜ_₁ * x + coeff_Θꜛₐ_ᵢₙₜ_₁ * _coeff_Θꜛₐ_ᵢₙₜ_2(T) * y - z 
                    end,
                    ∂Θₛₒₗ_∂Yₛₒₗ_int[1:M],∂ξₛₒₗ_∂Yₛₒₗ_int[1:M],T_int[1:M],∂Θₐ_∂Yₐ_int[1:M]) 

    half_len = Int(length(loss1) / 2) +1 
    for i in 1:length(loss1)
        λ = (i <= half_len) * i * (1.0 / half_len) + (i > half_len) * (2.0 - i / half_len)
        loss1[i] = λ * loss1[i]
    end

    loss = sum(abs2,loss1) / M  + sum(abs2,loss2) / M
    return loss
end
# --------------------------------------------------------------------------------------------
begin "test"
    ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1);  # allocation of ξₛₒₗ
    Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1);  # allocation of Θₛₒₗ
    Θꜛₐ = zeros(Float64,N + 2,M + 1);  # allocation of Θₐ
    ωꜛₐ = zeros(Float64,N + 2,M + 1);  # allocation of ωₐ
    ξ_int = zeros(Float64,M);
    T_int = zeros(Float64,M);
    Θꜛₐ_ᵢₙₜ = zeros(Float64,M);
    ωꜛ_int = zeros(Float64,M);
    B = zeros(Float64,N);
    α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N];
    α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N];

    Y = 0:Δy:1.0;
    @. Θꜛₛₒₗ[:,1] = -1.5 * Y^2 + 3.0 * Y ; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[:,1] .= 1.0;
    Θꜛₐ[:,1] .= 1.0; # Θₐ at x = 0
    ξꜛₛₒₗ[:,1] .= 1.0; # ξₛₒₗ at x = 0
    ωꜛₐ[:,1] .= 1.0; # ωₐ at x = 0
    
    X = vcat(0.99 * sort(rand(M);rev=true),sort(rand(M);rev=false));
    # X = vcat(-0.01 * sort(rand(M);rev=true),sort(rand(M);rev=true));

    # @code_warntype solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int)
    @time solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int,B,α_yₛₒₗ,α_yₐ)
end
X
heatmap(ξꜛₛₒₗ , title = "ξₛₒₗ - test")
heatmap(Θꜛₛₒₗ , title = "Θₛₒₗ - test")
heatmap(Θꜛₐ , title = "Θₐ - test")
heatmap(ωꜛₐ , title = "ωₐ - test")
# --------------------------------------------------------------------------------------------

#= begin "ForwardDiff gradient"
    der = zeros(Float64,2M);
    closure_dual = x -> ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}}(x, der...);

    function to_dual(x,f,der)
        x = ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64}}(x, der...)
        return x
    end

    function ∇f_f(X)
        ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of ξₛₒₗ
        Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of Θₛₒₗ
        Θꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of Θₐ
        ωꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of ωₐ
        ξ_int = zeros(Float64,M)
        # T_int = zeros(Float64,M)
        # Θꜛₐ_ᵢₙₜ = zeros(Float64,M)
        # ωꜛ_int = zeros(Float64,M)
        B = zeros(Float64,N)
        α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N]
        α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N]

        Θꜛₛₒₗ[:,1] .= 1.0 # Θₛₒₗ at x = 0
        Θꜛₐ[:,1] .= 1.0 # Θₐ at x = 0
        ξꜛₛₒₗ[:,1] .= 1.0 # ξₛₒₗ at x = 0
        ωꜛₐ[:,1] .= 1.0 # ωₐ at x = 0

        ξꜛₛₒₗ_Dual = closure_dual.(ξꜛₛₒₗ);
        ξ_int_Dual = closure_dual.(ξ_int);
        α_yₐ_Dual = closure_dual.(α_yₐ);
        α_yₛₒₗ_Dual = closure_dual.(α_yₛₒₗ);
        B_Dual = closure_dual.(B);

        Θꜛₛₒₗ_Dual = copy(ξꜛₛₒₗ_Dual);
        Θꜛₐ_Dual = copy(ξꜛₛₒₗ_Dual);
        ωꜛₐ_Dual = copy(ξꜛₛₒₗ_Dual);
        T_int_Dual = copy(ξ_int_Dual);
        Θꜛₐ_ᵢₙₜ_Dual = copy(ξ_int_Dual);
        ωꜛ_int_Dual = copy(ξ_int_Dual);
        der_ = Array(I(2M));
        der_X = [der_[:,i] for i in 1:2M];
        X_Dual = to_dual.(X,solve_domain!,der_X);
        sol_results = solve_domain!(X_Dual,ξꜛₛₒₗ_Dual,Θꜛₛₒₗ_Dual,Θꜛₐ_Dual,ωꜛₐ_Dual,
        ξ_int_Dual,T_int_Dual,Θꜛₐ_ᵢₙₜ_Dual,ωꜛ_int_Dual,B_Dual,α_yₛₒₗ_Dual,α_yₐ_Dual);

        f = sol_results.value
        ∇ = sol_results.partials
        return f,∇
    end

    X0 = vcat(0.9 * ones(Float64,M),sort(rand(M);rev=true))

    # @time f,∇ = ∇f_f(X);

    function f_(X)
        f,_ = ∇f_f(X)
        return f
    end

    function grad(X)
        _,∇ = ∇f_f(X)
        return ∇
    end
    # ∇ = zeros(Float64,2M)
    # f_(X0)
    # grad(X0)


    using Optim
    @time res = optimize(f_, grad, X0, BFGS();inplace = false)
    Optim.minimizer(res)[1:20]
    Optim.minimizer(res)[21:40]

    t = Δy:Δy:1.0
    u = Optim.minimizer(res)[M+1:2M]
    u = Optim.minimizer(res)[1:M]


    A = BSplineApprox(u, t, 2, 8, :ArcLen, :Average)
    scatter(t, u, label = "input data")
    plot!(A)
end =#
# --------------------------------------------------------------------------------------------

function _segment_solve!(X)
    ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of ξₛₒₗ
    Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of Θₛₒₗ
    Θꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of Θₐ
    ωꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of ωₐ
    ξ_int = zeros(Float64,M)
    T_int = zeros(Float64,M)
    Θꜛₐ_ᵢₙₜ = zeros(Float64,M)
    ωꜛ_int = zeros(Float64,M)

    B = zeros(Float64,N)
    α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N]
    α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N]

    Y = 0:Δy:1.0
    @. Θꜛₛₒₗ[:,1] = -1.5 * Y^2 + 3.0 * Y ; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[:,1] .= 1.0; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[1,1] = 0.0
    Θꜛₐ[:,1] .= 1.0 # Θₐ at x = 0
    ξꜛₛₒₗ[:,1] .= 1.0 # ξₛₒₗ at x = 0
    ωꜛₐ[:,1] .= 1.0 # ωₐ at x = 0

    loss = solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int,B,α_yₛₒₗ,α_yₐ)
    loss
end

X0 = vcat(1.1 * sort(rand(M);rev=true),1.1*sort(rand(M);rev=false));
@time _segment_solve!(X)
_segment_solve!(X0)
# options = Optim.Options(iterations= 1000,show_trace = true, show_every = 100)
# @time res = optimize(_segment_solve! ,X0, GradientDescent(),options;inplace = false)
# X0 = Optim.minimizer(res)
options = Optim.Options(iterations= 2000,show_trace = true, show_every = 100)
@time res = optimize(_segment_solve! ,X0, BFGS(),options;inplace = false)
# @time res = optimize(_segment_solve! ,X0, NewtonTrustRegion();inplace = false)
# X0 = Optim.minimizer(res)
# options = Optim.Options(iterations= 1000,show_trace = true, show_every = 100)
# @time res = optimize(_segment_solve! ,X0, NewtonTrustRegion(),options;inplace = false)
# X0 = Optim.minimizer(res)

u1 = Optim.minimizer(res)[1:M]
plot(u1, title = "ξₛₒₗ interface")

u2 = Optim.minimizer(res)[M+1:2M]
plot(u2, title = "Θₛₒₗ interface")


begin "test"
    ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1);  # allocation of ξₛₒₗ
    Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1);  # allocation of Θₛₒₗ
    Θꜛₐ = zeros(Float64,N + 2,M + 1);  # allocation of Θₐ
    ωꜛₐ = zeros(Float64,N + 2,M + 1);  # allocation of ωₐ
    ξ_int = zeros(Float64,M);
    T_int = zeros(Float64,M);
    Θꜛₐ_ᵢₙₜ = zeros(Float64,M);
    ωꜛ_int = zeros(Float64,M);
    B = zeros(Float64,N);
    α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N];
    α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N];

    Y = 0:Δy:1.0
    @. Θꜛₛₒₗ[:,1] = -1.5 * Y^2 + 3.0 * Y ; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[:,1] .= 1.0; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[1,1] = 0.0; # Θₛₒₗ at x = 0
    Θꜛₐ[:,1] .= 1.0; # Θₐ at x = 0
    ξꜛₛₒₗ[:,1] .= 1.0; # ξₛₒₗ at x = 0
    ωꜛₐ[:,1] .= 1.0; # ωₐ at x = 0
end
X = vcat(u1,u2)
# X = vcat(1.2 * ones(M),1.3 * ones(M))
# X[1:M] .= -X[1:M]
_segment_solve!(X)
# @code_warntype solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int)
@time solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int,B,α_yₛₒₗ,α_yₐ);

heatmap(ξꜛₛₒₗ , title = "ξₛₒₗ")
heatmap(Θꜛₛₒₗ , title = "Θₛₒₗ")
heatmap(Θꜛₐ , title = "Θₐ")
heatmap(ωꜛₐ , title = "ωₐ")

T_sol = zeros(Float64,N+2)
@. T_sol = Θꜛₛₒₗ[:,end] * (Tₛₒₗ_ᵢₙ - T_w) + T_w - 273.15
u_sol = zeros(Float64,N+2)
u_sol = Uꜛₛₒₗ.(collect(0:Δy:1.0)) .* Uₛₒₗ_ᵣ
T_sol_m = sum(T_sol .* u_sol) / sum(u_sol)

u_a = Uꜛₐ.(collect(0:Δy:1.0)) .* Uₐ_ᵣ
T_a = zeros(Float64,N+2)
@. T_a = Θꜛₐ[:,end] * (Tₐ_ᵢₙ - T_w) + T_w - 273.15
T_a_m = sum(T_a .* u_a) / sum(u_a)



mean(ξꜛₛₒₗ[:,end]) * ξₛₒₗ_ᵢₙ
mean(Θꜛₛₒₗ[:,end]) * (Tₛₒₗ_ᵢₙ - T_w) + T_w - 273.15
mean(Θꜛₐ[:,end]) * (Tₐ_ᵢₙ - T_w) + T_w - 273.15
(Θꜛₐ[:,end][end]) * (Tₐ_ᵢₙ - T_w) + T_w - 273.15
mean(ωꜛₐ[:,end]) * ωₐ_ᵢₙ
(ωꜛₐ[:,end][end]) * ωₐ_ᵢₙ


plot(Θꜛₛₒₗ[end ,:], title = "Θₛₒₗ interface")
plot(ξꜛₛₒₗ[end ,:], title = "ξꜛₛₒₗ interface")
plot(T_int .-273.15, title = "T_int")
plot(1 .- ξ_int, title = "ξ_int")
plot(ωꜛₐ[end,:] * ωₐ_ᵢₙ, title = "ωₐ")
ω_int = ωꜛₐ[end,1:M] .* ωₐ_ᵢₙ
1 .- ξ_int
T_int_r = reverse(T_int) .+ 6.0
plot(T_int_r .-273.15, title = "T_int_r")
plot(Θꜛₐ[:,end], title = "Θꜛₐ_out")
plot(ωꜛₐ[:,end], title = "ωₐ_out")

CoolProp.HAPropsSI("RH","T",17.73 + 273.15,"P",101325.0,"W",0.003655) * 100
first_der_at_int = (x , x⁻ , x⁻⁻) -> (3x - 4x⁻ + x⁻⁻) / 2Δy
∂ξₛₒₗ_∂Yₛₒₗ_int = map(first_der_at_int,ξꜛₛₒₗ[N+2,2:M+1],ξꜛₛₒₗ[N+1,2:M+1],ξꜛₛₒₗ[N,2:M+1])
∂Θₛₒₗ_∂Yₛₒₗ_int = map(first_der_at_int,Θꜛₛₒₗ[N+2,2:M+1],Θꜛₛₒₗ[N+1,2:M+1],Θꜛₛₒₗ[N,2:M+1])
∂Θₐ_∂Yₐ_int = map(first_der_at_int,Θꜛₐ[N+2,2:M+1],Θꜛₐ[N+1,2:M+1],Θꜛₐ[N,2:M+1]) |> reverse
∂ωₐ_∂Yₐ_int = map(first_der_at_int,ωꜛₐ[N+2,2:M+1],ωꜛₐ[N+1,2:M+1],ωꜛₐ[N,2:M+1]) |> reverse
loss1 = map((x,y) -> x + coeff_ωꜛₐ_ᵢₙₜ * y, ∂ωₐ_∂Yₐ_int[1:M],∂ξₛₒₗ_∂Yₛₒₗ_int[1:M])
loss2 = map((x,y,T,z) -> begin
                    - coeff_Θꜛₐ_ᵢₙₜ_₁ * x + coeff_Θꜛₐ_ᵢₙₜ_₁ * _coeff_Θꜛₐ_ᵢₙₜ_2(T) * y - z 
                end,
                ∂Θₛₒₗ_∂Yₛₒₗ_int[1:M],∂ξₛₒₗ_∂Yₛₒₗ_int[1:M],T_int[1:M],∂Θₐ_∂Yₐ_int[1:M]) 
plot(loss1, title = "loss1")
plot(loss2, title = "loss2")
plot(∂ξₛₒₗ_∂Yₛₒₗ_int, title = "∂ξₛₒₗ_∂Yₛₒₗ_int")
plot(∂Θₛₒₗ_∂Yₛₒₗ_int, title = "∂Θₛₒₗ_∂Yₛₒₗ_int")
plot(∂Θₐ_∂Yₐ_int, title = "∂Θₐ_∂Yₐ_int")
plot(∂ωₐ_∂Yₐ_int, title = "∂ωₐ_∂Yₐ_int")
plot(T_int .-273.15, title = "T_int")

plot(Θꜛₛₒₗ[:,1], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[:,2], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[:,3], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[:,4], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[:,5], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[:,6], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[:,7], title = "Θₛₒₗ -6")


plot(Θꜛₛₒₗ[:,end], title = "Θₛₒₗ")
plot(Θꜛₛₒₗ[:,end - 1], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[:,end - 2], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[:,end - 3], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[:,end - 4], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[:,end - 50], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[:,end - 60], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[:,end - 70], title = "Θₛₒₗ -6")

plot(Θꜛₛₒₗ[end ,:], title = "Θₛₒₗ0")
plot!(Θꜛₛₒₗ[end - 1,:], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[end - 2,:], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[end - 3,:], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[end - 4,:], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[end - 5,:], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[end - 6,:], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[end - 7,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[end - 8,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[end - 9,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[end - 10,:], title = "Θₛₒₗ -6")

plot(Θꜛₛₒₗ[  1,:], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[  2,:], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[  3,:], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[  4,:], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[  5,:], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[  6,:], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[  7,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[  8,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[  9,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[  10,:], title = "Θₛₒₗ -6")


plot(∂Θₛₒₗ_∂Yₛₒₗ_int[:], title = "∂Θₛₒₗ_∂Yₛₒₗ_int")
plot!(∂Θₐ_∂Yₐ_int[:], title = "∂Θₐ_∂Yₐ_int")

plot(Θꜛₛₒₗ[end ,:], title = "Θₛₒₗ interface")
plot(ξꜛₛₒₗ[end ,:], title = "ξꜛₛₒₗ interface")
plot(T_int .-273.15, title = "T_int")
plot(1 .- ξ_int, title = "ξ_int")
plot(ωꜛₐ[end,:] * ωₐ_ᵢₙ, title = "ωₐ")
ω_int = ωꜛₐ[end,1:M] .* ωₐ_ᵢₙ
1 .- ξ_int
T_int_r = reverse(T_int) .+ 6.0

contour(Θꜛₛₒₗ , title = "Θₛₒₗ",levels=1000, color=:turbo, clabels=false, cbar=true, lw=1)
contour(ξꜛₛₒₗ , title = "ξₛₒₗ",levels=1000, color=:turbo, clabels=false, cbar=true, lw=1)
contour(Θꜛₐ , title = "Θₐ",levels=1000, color=:turbo, clabels=false, cbar=true, lw=1)
contour(ωꜛₐ , title = "ωₐ",levels=1000, color=:turbo, clabels=false, cbar=true, lw=1)

i_eq(T,W) = (1.01 * (T -273.15) + W * (2500 + 1.84 * (T - 273.15))) * 1000
i_a(T,W) = CoolProp.HAPropsSI("H","T",T, "P",101325.0,"W",W)
plot(i_a.(T_int_r,ω_int) - i_eq.(T_int_r,ω_int))





p_v = _Pᵥₐₚₒᵣ_ₛₒₗ.(T_int,1.0 .- ξ_int);
plot(p_v[:], title = "p_v")
# --------------------------------------------------------------------------------------------

function _segment_solve!(W_B_NN)
    ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of ξₛₒₗ
    Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of Θₛₒₗ
    Θꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of Θₐ
    ωꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of ωₐ
    ξ_int = zeros(Float64,M)
    T_int = zeros(Float64,M)
    Θꜛₐ_ᵢₙₜ = zeros(Float64,M)
    ωꜛ_int = zeros(Float64,M)

    B = zeros(Float64,N)
    α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N]
    α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N]

    Y = 0:Δy:1.0
    @. Θꜛₛₒₗ[:,1] = -1.5 * Y^2 + 3.0 * Y ; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[:,1] .= 1.0; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[1,1] = 0.0
    Θꜛₐ[:,1] .= 1.0 # Θₐ at x = 0
    ξꜛₛₒₗ[:,1] .= 1.0 # ξₛₒₗ at x = 0
    ωꜛₐ[:,1] .= 1.0 # ωₐ at x = 0

    half_W_B_NN = Int(0.5 * length(W_B_NN))
    @unpack L1, L2, L3 = make_X_component_array(view(W_B_NN,1:half_W_B_NN),N_nr)
    simpleNN1(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)
    X1 = simpleNN1.(collect(Δy:Δy:1))

    @unpack L1, L2, L3 = make_X_component_array(view(W_B_NN,half_W_B_NN+1:length(W_B_NN)),N_nr)
    simpleNN2(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)
    X2 = simpleNN2.(collect(Δy:Δy:1))

    X = vcat(X1,X2)
    loss = solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int,B,α_yₛₒₗ,α_yₐ)
    # loss = loss + 1e6 * abs(simpleNN1(0.0) - 1.0)^2 + 1000.0 * abs(simpleNN2(0.0) - 1.5)^2
    loss
end

σ(x) = 1.0 / (1.0 + exp(-x))

function make_X_component_array(X,N)
    L1_ = (W = reshape(view(X,1:N),N,1), b = reshape(view(X,N+1:2N),N,1))
    L2_ = (W = reshape(view(X,2N+1:2N+N^2),N,N), b = reshape(view(X,2N+1+N^2:3N+N^2),N,1))
    L3_ = (W = reshape(view(X,3N+1+N^2:4 * N + N^2),1,N), b = reshape(view(X,4 * N + N^2 + 1),1))
    Layers = ComponentArray(L1 = L1_, L2 = L2_, L3 = L3_)
    return Layers
end


const N_nr = 2
ps_number = (4 * N_nr + N_nr^2 + 1) * 2

W_B_NN0 = rand(ps_number)
@time _segment_solve!(W_B_NN0)
# options = Optim.Options(x_tol = 1e-6, f_tol = 1e-3, g_tol = 1e-3, show_trace = true, show_every = 1)
# options = Optim.Options(iterations= 1000,show_trace = true, show_every = 100)
# @time res = optimize(_segment_solve! ,W_B_NN0, GradientDescent(),options;inplace = false)
# W_B_NN0 = Optim.minimizer(res)
options = Optim.Options(iterations= 2000,show_trace = true, show_every = 100)
@time res = optimize(_segment_solve! ,W_B_NN0, BFGS(),options;inplace = false)
W_B_NN0 = Optim.minimizer(res)
options = Optim.Options(iterations= 400,show_trace = true, show_every = 100)
@time res = optimize(_segment_solve! ,W_B_NN0, NewtonTrustRegion(),options;inplace = false)


u1_wb = Optim.minimizer(res)[1:Int(ps_number/2)]
@unpack L1, L2, L3 = make_X_component_array(u1_wb,N_nr)
simpleNN1(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)
u1 = simpleNN1.(collect(0:Δy:1.0))
plot(u1, title = "ξₛₒₗ interface")
simpleNN1(0.0)

u2_wb = Optim.minimizer(res)[Int(ps_number/2) + 1:ps_number]
@unpack L1, L2, L3 = make_X_component_array(u2_wb,N_nr)
simpleNN2(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)
u2 = simpleNN1.(collect(0:Δy:1.0))
# plot(u1, title = "ξₛₒₗ interface")
# u2 = Optim.minimizer(res)[M+1:2M]
plot(u2, title = "Θₛₒₗ interface")
simpleNN2(0.0)

begin "test"
    ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1);  # allocation of ξₛₒₗ
    Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1);  # allocation of Θₛₒₗ
    Θꜛₐ = zeros(Float64,N + 2,M + 1);  # allocation of Θₐ
    ωꜛₐ = zeros(Float64,N + 2,M + 1);  # allocation of ωₐ
    ξ_int = zeros(Float64,M);
    T_int = zeros(Float64,M);
    Θꜛₐ_ᵢₙₜ = zeros(Float64,M);
    ωꜛ_int = zeros(Float64,M);
    B = zeros(Float64,N);
    α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N];
    α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N];

    Y = 0:Δy:1.0
    @. Θꜛₛₒₗ[:,1] = -1.5 * Y^2 + 3.0 * Y ; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[:,1] .= 1.0; # Θₛₒₗ at x = 0
    # Θꜛₛₒₗ[1,1] = 0.0; # Θₛₒₗ at x = 0
    Θꜛₐ[:,1] .= 1.0; # Θₐ at x = 0
    ξꜛₛₒₗ[:,1] .= 1.0; # ξₛₒₗ at x = 0
    ωꜛₐ[:,1] .= 1.0; # ωₐ at x = 0
end
X = vcat(u1,u2)
# X[1:M] .= -X[1:M]
_segment_solve!(X)
# @code_warntype solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int)
@time solve_domain!(X,ξꜛₛₒₗ,Θꜛₛₒₗ,Θꜛₐ,ωꜛₐ,ξ_int,T_int,Θꜛₐ_ᵢₙₜ,ωꜛ_int,B,α_yₛₒₗ,α_yₐ);

heatmap(ξꜛₛₒₗ , title = "ξₛₒₗ")
heatmap(Θꜛₛₒₗ , title = "Θₛₒₗ")
heatmap(Θꜛₐ , title = "Θₐ")
heatmap(ωꜛₐ , title = "ωₐ")

T_sol = zeros(Float64,N+2)
@. T_sol = Θꜛₛₒₗ[:,end] * (Tₛₒₗ_ᵢₙ - T_w) + T_w - 273.15
u_sol = zeros(Float64,N+2)
u_sol = Uꜛₐ.(collect(0:Δy:1.0)) .* Uₛₒₗ_ᵣ
T_sol_m = sum(T_sol .* u_sol) / sum(u_sol)

mean(ξꜛₛₒₗ[:,end]) * ξₛₒₗ_ᵢₙ
mean(Θꜛₛₒₗ[:,end]) * (Tₛₒₗ_ᵢₙ - T_w) + T_w - 273.15
mean(Θꜛₐ[:,end]) * (Tₐ_ᵢₙ - T_w) + T_w - 273.15
mean(ωꜛₐ[:,end]) * ωₐ_ᵢₙ
first_der_at_int = (x , x⁻ , x⁻⁻) -> (3x - 4x⁻ + x⁻⁻) / 2Δy
∂ξₛₒₗ_∂Yₛₒₗ_int = map(first_der_at_int,ξꜛₛₒₗ[N+2,2:M+1],ξꜛₛₒₗ[N+1,2:M+1],ξꜛₛₒₗ[N,2:M+1])
∂Θₛₒₗ_∂Yₛₒₗ_int = map(first_der_at_int,Θꜛₛₒₗ[N+2,2:M+1],Θꜛₛₒₗ[N+1,2:M+1],Θꜛₛₒₗ[N,2:M+1])
∂Θₐ_∂Yₐ_int = map(first_der_at_int,Θꜛₐ[N+2,2:M+1],Θꜛₐ[N+1,2:M+1],Θꜛₐ[N,2:M+1]) |> reverse
∂ωₐ_∂Yₐ_int = map(first_der_at_int,ωꜛₐ[N+2,2:M+1],ωꜛₐ[N+1,2:M+1],ωꜛₐ[N,2:M+1]) |> reverse
loss1 = map((x,y) -> x + coeff_ωꜛₐ_ᵢₙₜ * y, ∂ωₐ_∂Yₐ_int[1:M],∂ξₛₒₗ_∂Yₛₒₗ_int[1:M])
loss2 = map((x,y,T,z) -> begin
                    - coeff_Θꜛₐ_ᵢₙₜ_₁ * x + coeff_Θꜛₐ_ᵢₙₜ_₁ * _coeff_Θꜛₐ_ᵢₙₜ_2(T) * y - z 
                end,
                ∂Θₛₒₗ_∂Yₛₒₗ_int[1:M],∂ξₛₒₗ_∂Yₛₒₗ_int[1:M],T_int[1:M],∂Θₐ_∂Yₐ_int[1:M]) 
plot(loss1, title = "loss1")
plot(loss2, title = "loss2")
plot(∂ξₛₒₗ_∂Yₛₒₗ_int, title = "∂ξₛₒₗ_∂Yₛₒₗ_int")
plot(∂Θₛₒₗ_∂Yₛₒₗ_int, title = "∂Θₛₒₗ_∂Yₛₒₗ_int")
plot(∂Θₐ_∂Yₐ_int, title = "∂Θₐ_∂Yₐ_int")
plot(∂ωₐ_∂Yₐ_int, title = "∂ωₐ_∂Yₐ_int")
plot(T_int .-273.15, title = "T_int")

plot(Θꜛₛₒₗ[:,1], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[:,2], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[:,3], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[:,4], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[:,5], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[:,6], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[:,7], title = "Θₛₒₗ -6")

plot(Θꜛₛₒₗ[end ,:], title = "Θₛₒₗ0")
plot!(Θꜛₛₒₗ[end - 1,:], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[end - 2,:], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[end - 3,:], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[end - 4,:], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[end - 5,:], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[end - 6,:], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[end - 7,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[end - 8,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[end - 9,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[end - 10,:], title = "Θₛₒₗ -6")

plot(Θꜛₛₒₗ[  1,:], title = "Θₛₒₗ")
plot!(Θꜛₛₒₗ[  2,:], title = "Θₛₒₗ -1")
plot!(Θꜛₛₒₗ[  3,:], title = "Θₛₒₗ -2")
plot!(Θꜛₛₒₗ[  4,:], title = "Θₛₒₗ -3")
plot!(Θꜛₛₒₗ[  5,:], title = "Θₛₒₗ -4")
plot!(Θꜛₛₒₗ[  6,:], title = "Θₛₒₗ -5")
plot!(Θꜛₛₒₗ[  7,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[  8,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[  9,:], title = "Θₛₒₗ -6")
plot!(Θꜛₛₒₗ[  10,:], title = "Θₛₒₗ -6")


plot(∂Θₛₒₗ_∂Yₛₒₗ_int[20:50], title = "∂Θₛₒₗ_∂Yₛₒₗ_int")
plot!(∂Θₐ_∂Yₐ_int[20:50], title = "∂Θₐ_∂Yₐ_int")

plot(Θꜛₛₒₗ[end ,:], title = "Θₛₒₗ interface")
plot(ξꜛₛₒₗ[end ,:], title = "ξꜛₛₒₗ interface")
plot(T_int .-273.15, title = "T_int")
plot(1 .- ξ_int, title = "ξ_int")
p_v = _Pᵥₐₚₒᵣ_ₛₒₗ.(T_int,1.0 .- ξ_int)
plot(p_v[20:80], title = "p_v")
# --------------------------------------------------------------------------------------------

CoolProp.HAPropsSI("R","P",101325.0,"W",0.0096,"T",273.15 + 17.85)

# --------------------------------------------------------------------------------------------
using JLD2

@load "u1.jld2" u1
@load "u2.jld2" u2
reverse!(u1)
push!(u1,1.0)
reverse!(u1)
using Lux , Random, OptimizationOptimisers, Printf,ComponentArrays
using Parameters: @unpack
σ(x) = 1.0 / (1.0 + exp(-x))
N_nerouns = 3
function make_X_component_array(X,N)
    L1_ = (W = reshape(view(X,1:N),N,1), b = reshape(view(X,N+1:2N),N,1))
    L2_ = (W = reshape(view(X,2N+1:2N+N^2),N,N), b = reshape(view(X,2N+1+N^2:3N+N^2),N,1))
    L3_ = (W = reshape(view(X,3N+1+N^2:4 * N_nerouns + N_nerouns^2),1,N), b = reshape(view(X,4 * N_nerouns + N_nerouns^2 + 1),1))
    Layers = ComponentArray(L1 = L1_, L2 = L2_, L3 = L3_)
    return Layers
end
X = rand(4 * N_nerouns + N_nerouns^2 + 1)
@unpack L1, L2, L3 = make_X_component_array(X,N_nerouns)
simpleNN(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)

function loss(X,y)
    @unpack L1, L2, L3 = make_X_component_array(X,N_nerouns)
    simpleNN(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)
    ŷ = simpleNN.(collect(0:Δy:1))
    return sum(abs2,y - ŷ)
end

loss_closure = ps -> loss(ps,u1)
X0_ = rand(4 * N_nerouns + N_nerouns^2 + 1)
options_ = Optim.Options(f_abstol=1e-6,g_abstol=1e-6,iterations=2_000,show_trace = true, show_every = 100)
res_model = optimize(loss_closure,X0_,BFGS(),options_;inplace = false)
X_model = Optim.minimizer(res_model)
@unpack L1, L2, L3 = make_X_component_array(X_model,N_nerouns)
simpleNN(x) = first(L3.W*σ.(L2.W*σ.(L1.W*x + L1.b) + L2.b) + L3.b)
plot(simpleNN.(collect(Δy:Δy:1)))
plot!(u1)

# --------------------------------------------------------------------------------------------

begin
    ξꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of ξₛₒₗ
    Θꜛₛₒₗ = zeros(Float64,N + 2,M + 1)  # allocation of Θₛₒₗ
    Θꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of Θₐ
    ωꜛₐ = zeros(Float64,N + 2,M + 1)  # allocation of ωₐ
    ξ_int = zeros(Float64,M)
    T_int = zeros(Float64,M)
    Θꜛₐ_ᵢₙₜ = zeros(Float64,M)
    ωꜛ_int = zeros(Float64,M)
    B = zeros(Float64,N)
    α_yₐ = [Uꜛₐ(i * Δy) for i in 1:N]
    α_yₛₒₗ = [Uꜛₛₒₗ(i * Δy) for i in 1:N]

    Θꜛₛₒₗ[:,1] .= 1.0 # Θₛₒₗ at x = 0
    Θꜛₐ[:,1] .= 1.0 # Θₐ at x = 0
    ξꜛₛₒₗ[:,1] .= 1.0 # ξₛₒₗ at x = 0
    ωꜛₐ[:,1] .= 1.0 # ωₐ at x = 0
end

begin "X initiating"
    X = vcat(0.9 * ones(Float64,M),sort(rand(M);rev=true))
end

begin "ForwardDiff Dulaize the variables initiating"
    der = zeros(Float64,2M);
    closure_dual = x -> ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}}(x, der...);
    ξꜛₛₒₗ_Dual_ = closure_dual.(ξꜛₛₒₗ);
    ξ_int_Dual_ = closure_dual.(ξ_int);
    α_yₐ_Dual = closure_dual.(α_yₐ);
    α_yₛₒₗ_Dual = closure_dual.(α_yₛₒₗ);
    B_Dual_ = closure_dual.(B);
end

begin "copy dual variables"
    ξꜛₛₒₗ_Dual = copy(ξꜛₛₒₗ_Dual_);
    Θꜛₛₒₗ_Dual = copy(ξꜛₛₒₗ_Dual);
    Θꜛₐ_Dual = copy(ξꜛₛₒₗ_Dual);
    ωꜛₐ_Dual = copy(ξꜛₛₒₗ_Dual);
    ξ_int_Dual = copy(ξ_int_Dual_);
    T_int_Dual = copy(ξ_int_Dual);
    Θꜛₐ_ᵢₙₜ_Dual = copy(ξ_int_Dual);
    ωꜛ_int_Dual = copy(ξ_int_Dual);
    B_Dual = copy(B_Dual_);
end

begin "Dualize the X i.e. guess ξ and Θ @ interface nodes"
    der_ = Array(I(2M));
    der_X = [der_[:,i] for i in 1:2M];
    function to_dual(x,f,der)
        x = ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64}}(x, der...)
        return x
    end
    X_Dual = to_dual.(X,solve_domain!,der_X);
end

@time sol_results = solve_domain!(X_Dual,ξꜛₛₒₗ_Dual,Θꜛₛₒₗ_Dual,Θꜛₐ_Dual,ωꜛₐ_Dual,
    ξ_int_Dual,T_int_Dual,Θꜛₐ_ᵢₙₜ_Dual,ωꜛ_int_Dual,B_Dual,α_yₛₒₗ_Dual,α_yₐ_Dual);

loss = sol_results.value
∇ = sol_results.partials[1:2M]

# --------------------------------------------------------------------------------------------       
begin "hessian"
    gate = false
    if gate
        make_dual =  ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}} # make dual variable
        der = zeros(Float64,2M);
        zero_tuple_dual = make_dual(0.0, der...);
        zero_tuple_duals = [zero_tuple_dual for i in 1:2M+1];
        closure_dual_hessian = ξꜛₛₒₗ_val -> make_dual(make_dual(ξꜛₛₒₗ_val, der...), zero_tuple_duals...)

        begin "ForwardDiff Dulaize the variables initiating"
            ξꜛₛₒₗ_Dual_ = closure_dual_hessian.(ξꜛₛₒₗ);
            ξ_int_Dual_ = closure_dual_hessian.(ξ_int);
            α_yₐ_Dual = closure_dual_hessian.(α_yₐ);
            α_yₛₒₗ_Dual = closure_dual_hessian.(α_yₛₒₗ);
            B_Dual_ = closure_dual_hessian.(B);
        end 

        ξꜛₛₒₗ_Dual = copy(ξꜛₛₒₗ_Dual_);
        Θꜛₛₒₗ_Dual = copy(ξꜛₛₒₗ_Dual);
        Θꜛₐ_Dual = copy(ξꜛₛₒₗ_Dual);
        ωꜛₐ_Dual = copy(ξꜛₛₒₗ_Dual);
        ξ_int_Dual = copy(ξ_int_Dual_);
        T_int_Dual = copy(ξ_int_Dual);
        Θꜛₐ_ᵢₙₜ_Dual = copy(ξ_int_Dual);
        ωꜛ_int_Dual = copy(ξ_int_Dual);
        B_Dual = copy(B_Dual_);

        der = Array(I(2M));
        X_dual = Vector{ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}, ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}, Float64, 2M}, 2M}}(undef,2M);

        pprint_struct(X_dual[2].value.partials)

        x_grad_dual = [make_dual(X[i], der[:,i]...) for i in 1:2M];
        unit_tuple_dual = make_dual(1.0, zeros(Float64,2M)...);
        for i in 1:2M
            tmp = der[:,i] .* unit_tuple_dual
            X_dual[i] = make_dual(make_dual(X[i], der[:,i]...),tmp...)
        end

        @time sol_results = solve_domain!(X_dual,ξꜛₛₒₗ_Dual,Θꜛₛₒₗ_Dual,Θꜛₐ_Dual,ωꜛₐ_Dual,
            ξ_int_Dual,T_int_Dual,Θꜛₐ_ᵢₙₜ_Dual,ωꜛ_int_Dual,B_Dual,α_yₛₒₗ_Dual,α_yₐ_Dual);
    end
end
# --------------------------------------------------------------------------------------------
begin "hessian demo / MWE"
    using ForwardDiff
    using GarishPrint

    gcubic(x,y) = begin 

        x[1]*x[1]*x[2] + x[1] + x[2] * y[2] + x[3] * y[1]
    end
    # xx = ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}}(1, 1,0)
    # yy = ForwardDiff.Dual{ForwardDiff.Tag{typeof(solve_domain!), Float64}}(2, 0,1)


    x1    =  ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(3.0,1.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(1.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0));

    x2    = ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(2.0,0.0,1.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(1.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0));

    x3    = ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(4.0,0.0,0.0,1.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(1.0,0.0,0.0,0.0));
    x = [x1,x2,x3];
    y1   =  ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(10.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0));

    y2    = ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(20.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0),
            ForwardDiff.Dual{ForwardDiff.Tag{typeof(gcubic), Float64}}(0.0,0.0,0.0,0.0));

    y = [y1 y2 y1;
        y1 y2 y2;
        y1 y2 y1;]

    res = gcubic(x,y)
end



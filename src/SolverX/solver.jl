export solve_coil_ode
# export solve_coil_naive_ode

include("NTU.jl")

function ionic_liquid_coil_ode!(du,u, p, t)
    # ωₐᵢᵣ, iₐᵢᵣ, ṁₛₒₗ,ξₛₒₗ, iₛₒₗ = u
    # ========================================
    Le = p[1]
    ∂Qᵣ = p[2]
    ṁₐᵢᵣ = p[3]
    NTUᴰₐᵢᵣ = p[4]
    σ = p[5]
    ṁₛₒₗ_ᵢₙ =  p[6]
    ξₛₒₗ_ᵢₙ = p[7]
    iₛₒₗ_ᵢₙ = p[8]
    ωₐ_ᵢₙ = p[9]
    iₐ_ᵢₙ = p[10]
    MR = ṁₛₒₗ_ᵢₙ / ṁₐᵢᵣ
    ER = iₛₒₗ_ᵢₙ / iₐ_ᵢₙ
    # ========================================
    Tₛₒₗ = calculate_T_sol(u[5] * iₛₒₗ_ᵢₙ, u[4] * ξₛₒₗ_ᵢₙ, IL)
    Pᵥₐₚₒᵣ_ₛₒₗ = _Pᵥₐₚₒᵣ_ₛₒₗ(Tₛₒₗ,u[4] * ξₛₒₗ_ᵢₙ,IL)
    ωₑ = 0.622 * Pᵥₐₚₒᵣ_ₛₒₗ / (101325.0 - Pᵥₐₚₒᵣ_ₛₒₗ) / ωₐ_ᵢₙ
    iₑ = (1.01 * (Tₛₒₗ - 273.15) + ωₑ * ωₐ_ᵢₙ * (2500 + 1.04 * (Tₛₒₗ - 273.15))) / iₐ_ᵢₙ
    iₑ *= 1000
    iᵥₐₚₒᵣ_ₜₛ = iᵥ_ₛₐₜ(Tₛₒₗ) / iₐ_ᵢₙ 

    du[1] = σ * NTUᴰₐᵢᵣ * (u[1] - ωₑ)
    du[2] = σ * NTUᴰₐᵢᵣ * Le * ((u[2] - iₑ) + (ωₐ_ᵢₙ * iᵥₐₚₒᵣ_ₜₛ * (1 / Le - 1) * (u[1] - ωₑ)))
    du[3] = σ * ωₐ_ᵢₙ * du[1] / MR
    du[4] = (-u[4] / u[3]) * du[3]
    du[5] = (1 / u[3]) * (σ * (1.0 / (MR * ER)) * du[2] - u[5] * du[3] - ∂Qᵣ / (ṁₛₒₗ_ᵢₙ * iₛₒₗ_ᵢₙ))
    nothing
end

function bca!(res_a, u_a, p)
    res_a[1] = u_a[3] - 1.0
    res_a[2] = u_a[4] - 1.0
    res_a[3] = u_a[5] - 1.0
    nothing
end

function bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1.0
    res_b[2] = u_b[2] - 1.0
    nothing
end

function solve_coil_ode(IL, Le ,∂Qᵣ ,ṁₐᵢᵣ_ᵢₙ ,NTUᴰₐᵢᵣ ,σ ,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ,
                dt,tspan)
    IonicLiquidDeh1DSteady.IL = IL
    p = @SVector[Le, ∂Qᵣ, ṁₐᵢᵣ_ᵢₙ, NTUᴰₐᵢᵣ, σ, ṁₛₒₗ_ᵢₙ, ξₛₒₗ_ᵢₙ, iₛₒₗ_ᵢₙ, ωₐ_ᵢₙ, iₐ_ᵢₙ]
    u0 = [0.1, 0.1 , 1.0001 , 0.9 , 1.01]

    bvp_fun = BVPFunction(
                    ionic_liquid_coil_ode!, (bca!, bcb!);
                    bcresid_prototype = (zeros(3), zeros(2)), twopoint=Val(true)
        )

    prob = TwoPointBVProblem(bvp_fun,
                u0,
                tspan,
                p)
    sol = solve(prob, MIRK6(), dt = dt)

    len_vec = _len_sol(sol.u)
    
    ωₐᵢᵣ = Vector{Float64}(undef,len_vec)
    iₐᵢᵣ = Vector{Float64}(undef,len_vec)
    ṁₛₒₗ = Vector{Float64}(undef,len_vec)
    ξₛₒₗ = Vector{Float64}(undef,len_vec)
    iₛₒₗ = Vector{Float64}(undef,len_vec)
    t = Vector{Float64}(undef,len_vec)
    _assign_data!(sol.u,sol.t,t,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ,len_vec,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ)

    t,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ
end

_len_sol(x::Vector) = length(x)

function _assign_data!(data::Vector{Vector{T}},t_::Vector{T},t::Vector{T},ωₐᵢᵣ::Vector{T}, iₐᵢᵣ::Vector{T},
     ṁₛₒₗ::Vector{T}, ξₛₒₗ::Vector{T}, iₛₒₗ::Vector{T},len_vec,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ) where T<:AbstractFloat
    @inbounds for i in 1:len_vec
        # ωₐᵢᵣ, iₐᵢᵣ, ṁₛₒₗ,ξₛₒₗ, iₛₒₗ = u
        t[i] = t_[i]
        ωₐᵢᵣ[i] = data[i][1] * ωₐ_ᵢₙ
        iₐᵢᵣ[i] = data[i][2] * iₐ_ᵢₙ
        ṁₛₒₗ[i] = data[i][3] * ṁₛₒₗ_ᵢₙ
        ξₛₒₗ[i] = data[i][4] * ξₛₒₗ_ᵢₙ
        iₛₒₗ[i] = data[i][5] * iₛₒₗ_ᵢₙ
    end
    nothing
end


# ========================================
#= function solve_coil_ode!(IL ,H ,Le ,∂Qᵣ ,ṁₐᵢᵣ_ᵢₙ ,NTUᴰₐᵢᵣ ,σ ,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ,
                        dt,tspan,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ)
    p = @SVector[IL, H, Le, ∂Qᵣ, ṁₐᵢᵣ_ᵢₙ, NTUᴰₐᵢᵣ, σ, ṁₛₒₗ_ᵢₙ, ξₛₒₗ_ᵢₙ, iₛₒₗ_ᵢₙ, ωₐ_ᵢₙ, iₐ_ᵢₙ]
    u0 = [0.1, 0.1 , 1.0001 , 0.9 , 1.01]
    
    bvp_fun = BVPFunction(
                ionic_liquid_coil_ode!, (bca!, bcb!);
                bcresid_prototype = (zeros(3), zeros(2)), twopoint = Val(true)
                )
    
    prob = TwoPointBVProblem(bvp_fun,
        u0,
        tspan,
        p)
    
    sol = solve(prob, MIRK4(), dt = dt)
    # sol = solve(prob,
    #             MultipleShooting(
    #                             10, Rodas5(), NewtonRaphson(; autodiff = AutoForwardDiff(; chunksize = 2)))
    #                             , dt = dt)
    # sol = solve(prob, Shooting(Tsit5()))
        
    @inbounds for i in 1:length(sol.u)
        # ωₐᵢᵣ, iₐᵢᵣ, ṁₛₒₗ,ξₛₒₗ, iₛₒₗ = u
        ωₐᵢᵣ[i] = sol.u[i][1] * ωₐ_ᵢₙ
        iₐᵢᵣ[i] = sol.u[i][2] * iₐ_ᵢₙ
        ṁₛₒₗ[i] = sol.u[i][3] * ṁₛₒₗ_ᵢₙ
        ξₛₒₗ[i] = sol.u[i][4] * ξₛₒₗ_ᵢₙ
        iₛₒₗ[i] = sol.u[i][5] * iₛₒₗ_ᵢₙ
    end
    nothing
end =#

# ========================================
#= 
function ionic_liquid_coil_naive_ode!(du,u, p, t)
    # ωₐᵢᵣ, iₐᵢᵣ, ṁₛₒₗ,ξₛₒₗ, iₛₒₗ = u
    # ========================================
    IL = p[1]
    H = p[2]
    Le = p[3]
    ∂Qᵣ = p[4]
    ṁₐᵢᵣ = p[5]
    NTUᴰₐᵢᵣ = p[6]
    σ = p[7]
    ṁₛₒₗ_ᵢₙ =  p[8]
    ξₛₒₗ_ᵢₙ = p[9]
    iₛₒₗ_ᵢₙ = p[10]
    ωₐ_ᵢₙ = p[11]
    iₐ_ᵢₙ = p[12]
    MR = ṁₛₒₗ_ᵢₙ / ṁₐᵢᵣ
    ER = iₛₒₗ_ᵢₙ / iₐ_ᵢₙ
    # ========================================
    Tₛₒₗ = calculate_T_sol(u[3] * iₛₒₗ_ᵢₙ, ξₛₒₗ_ᵢₙ, IL)
    Pᵥₐₚₒᵣ_ₛₒₗ = _Pᵥₐₚₒᵣ_ₛₒₗ(Tₛₒₗ, ξₛₒₗ_ᵢₙ,IL)
    ωₑ = 0.622 * Pᵥₐₚₒᵣ_ₛₒₗ / (101325.0 - Pᵥₐₚₒᵣ_ₛₒₗ) / ωₐ_ᵢₙ
    iₑ = (1.01 * (Tₛₒₗ - 273.15) + ωₑ * ωₐ_ᵢₙ * (2500 + 1.04 * (Tₛₒₗ - 273.15))) / iₐ_ᵢₙ
    iₑ *= 1000
    iᵥₐₚₒᵣ_ₜₛ = iᵥ_ₛₐₜ(Tₛₒₗ) / iₐ_ᵢₙ 

    du[1] = σ * NTUᴰₐᵢᵣ * (u[1] - ωₑ)
    du[2] = σ * NTUᴰₐᵢᵣ * Le * ((u[2] - iₑ) + (ωₐ_ᵢₙ * iᵥₐₚₒᵣ_ₜₛ * (1 / Le - 1) * (u[1] - ωₑ)))
    du[3] = (1 / ṁₛₒₗ_ᵢₙ) * (σ * (1.0 / (MR * ER)) * du[2] - ∂Qᵣ / (ṁₛₒₗ_ᵢₙ * iₛₒₗ_ᵢₙ))
    nothing
end

function bca_naive!(res_a, u_a, p)
    res_a[1] = u_a[3] - 1.0
    nothing
end

function bcb_naive!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1.0
    res_b[2] = u_b[2] - 1.0
    nothing
end

function solve_coil_naive_ode(IL ,H ,Le ,∂Qᵣ ,ṁₐᵢᵣ_ᵢₙ ,NTUᴰₐᵢᵣ ,σ ,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ,
        dt,tspan)
    p = @SVector[IL, H, Le, ∂Qᵣ, ṁₐᵢᵣ_ᵢₙ, NTUᴰₐᵢᵣ, σ, ṁₛₒₗ_ᵢₙ, ξₛₒₗ_ᵢₙ, iₛₒₗ_ᵢₙ, ωₐ_ᵢₙ, iₐ_ᵢₙ]
    u0 = [0.1, 0.1  , 1.01]

    bvp_fun = BVPFunction(
        ionic_liquid_coil_naive_ode!, (bca_naive!, bcb_naive!);
        bcresid_prototype = (zeros(1), zeros(2)), twopoint = Val(true)
    )

    prob = TwoPointBVProblem(bvp_fun,
        u0,
        tspan,
        p)

    sol = solve(prob, MIRK4(defect_threshold = 0.1, max_num_subintervals = 30000), dt = dt)

    len_vec = length(sol.u)
    ωₐᵢᵣ = zeros(len_vec)
    iₐᵢᵣ = zeros(len_vec)
    iₛₒₗ = zeros(len_vec)
    @inbounds for i in 1:length(sol.u)
        # ωₐᵢᵣ, iₐᵢᵣ, iₛₒₗ = u
        ωₐᵢᵣ[i] = sol.u[i][1] * ωₐ_ᵢₙ
        iₐᵢᵣ[i] = sol.u[i][2] * iₐ_ᵢₙ
        iₛₒₗ[i] = sol.u[i][3] * iₛₒₗ_ᵢₙ
    end
    ωₐᵢᵣ,iₐᵢᵣ,iₛₒₗ
end =#
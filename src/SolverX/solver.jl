export solve_coil_ode


function ionic_liquid_coil_ode!(du,u, p, t)
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
    Tₛₒₗ = calculate_T_sol(u[5] * iₛₒₗ_ᵢₙ, u[4] * ξₛₒₗ_ᵢₙ, IL)
    # @info "Tₛₒₗ = $(Tₛₒₗ - 273.15)"
    Pᵥₐₚₒᵣ_ₛₒₗ = _Pᵥₐₚₒᵣ_ₛₒₗ(Tₛₒₗ,u[4] * ξₛₒₗ_ᵢₙ,IL)
    # @info "Pᵥₐₚₒᵣ_ₛₒₗ = $(Pᵥₐₚₒᵣ_ₛₒₗ)"
    ωₑ = 0.622 * Pᵥₐₚₒᵣ_ₛₒₗ / (101325.0 - Pᵥₐₚₒᵣ_ₛₒₗ) / ωₐ_ᵢₙ
    iₑ = (1.01 * (Tₛₒₗ - 273.15) + ωₑ * ωₐ_ᵢₙ * (2500 + 1.04 * (Tₛₒₗ - 273.15))) / iₐ_ᵢₙ
    iₑ *= 1000
    iᵥₐₚₒᵣ_ₜₛ = iᵥ_ₛₐₜ(Tₛₒₗ) / iₐ_ᵢₙ 
    # @info "iᵥₐₚₒᵣ_ₜₛ = $(iᵥₐₚₒᵣ_ₜₛ)"

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

function solve_coil_ode(IL ,H ,Le ,∂Qᵣ ,ṁₐᵢᵣ_ᵢₙ ,NTUᴰₐᵢᵣ ,σ ,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ)
    p = @SVector[IL, H, Le, ∂Qᵣ, ṁₐᵢᵣ_ᵢₙ, NTUᴰₐᵢᵣ, σ, ṁₛₒₗ_ᵢₙ, ξₛₒₗ_ᵢₙ, iₛₒₗ_ᵢₙ, ωₐ_ᵢₙ, iₐ_ᵢₙ]
    u0 = [0.7, 0.7 , 1.0001 , 0.9 , 1.01 ]
    bvp_fun = BVPFunction(
                ionic_liquid_coil_ode!, (bca!, bcb!),
                bcresid_prototype = (zeros(3), zeros(2)), twopoint = Val(true)
                )
    tspan = (0.0, 1.0)
    prob = BVProblem(bvp_fun,
        u0,
        tspan,
        p)
    sol = solve(prob, MIRK4(), dt = 0.001)
    sol
end






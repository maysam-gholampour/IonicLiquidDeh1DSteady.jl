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
    # ========================================
    Tₛₒₗ = calculate_T_sol(u[5],u[4], IL)
    Pᵥₐₚₒᵣ_ₛₒₗ = _Pᵥₐₚₒᵣ_ₛₒₗ(Tₛₒₗ,u[4],IL)
    ωₑ = 0.622 * Pᵥₐₚₒᵣ_ₛₒₗ / (101325.0 - Pᵥₐₚₒᵣ_ₛₒₗ)
    iₑ = 1.01 * (Tₛₒₗ - 273.15) + ωₑ * (2500 + 1.04 * (Tₛₒₗ - 273.15))
    iₑ *= 1000
    iᵥₐₚₒᵣ_ₜₛ = iᵥ_ₛₐₜ(Tₛₒₗ)

    du[1] = σ * (NTUᴰₐᵢᵣ / H) * (u[1] - ωₑ)
    du[2] = σ * (NTUᴰₐᵢᵣ * Le / H) * ((u[2] - iₑ) + (iᵥₐₚₒᵣ_ₜₛ * (1 / Le - 1) * (u[1] - ωₑ)))
    du[3] = σ * ṁₐᵢᵣ * du[1]
    du[4] = (-u[4] / u[3]) * du[3]
    du[5] = (1 / u[3]) * (σ * ṁₐᵢᵣ * du[2] - u[5] * du[3] - ∂Qᵣ / H)
    nothing
end

function bca!(res_a, u_a, p)
    ṁₛₒₗ_ᵢₙ =  p[8]
    ξₛₒₗ_ᵢₙ = p[9]
    iₛₒₗ_ᵢₙ = p[10]
    res_a[1] = u_a[3] - ṁₛₒₗ_ᵢₙ
    res_a[2] = u_a[4] - ξₛₒₗ_ᵢₙ
    res_a[3] = u_a[5] - iₛₒₗ_ᵢₙ
    nothing
end

function bcb!(res_b, u_b, p)
    ωₐ_ᵢₙ = p[11]
    iₐ_ᵢₙ = p[12]
    res_b[1] = u_b[1] - ωₐ_ᵢₙ
    res_b[2] = u_b[2] - iₐ_ᵢₙ
    nothing
end

function solve_coil_ode(IL ,H ,Le ,∂Qᵣ ,ṁₐᵢᵣ_ᵢₙ ,NTUᴰₐᵢᵣ ,σ ,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ)
    p = @SVector[IL, H, Le, ∂Qᵣ, ṁₐᵢᵣ_ᵢₙ, NTUᴰₐᵢᵣ, σ, ṁₛₒₗ_ᵢₙ, ξₛₒₗ_ᵢₙ, iₛₒₗ_ᵢₙ, ωₐ_ᵢₙ, iₐ_ᵢₙ]
    u0 = [0.7 * ωₐ_ᵢₙ, 0.7 * iₐ_ᵢₙ, 1.0001 * ṁₛₒₗ_ᵢₙ, 0.9 *ξₛₒₗ_ᵢₙ, 1.01 * iₛₒₗ_ᵢₙ]
    bvp_fun = BVPFunction(
                ionic_liquid_coil_ode!, (bca!, bcb!),
                bcresid_prototype = (zeros(3), zeros(2)), twopoint = Val(true)
                )
    tspan = (0.0, H)
    prob = BVProblem(bvp_fun,
        u0,
        tspan,
        p)
    sol = solve(prob, MIRK4(), dt = 0.001)
    sol
end






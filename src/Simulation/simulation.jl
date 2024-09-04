export simulate!

function simulate!(plateFinCircularTube::PlateFinCircularTube,fluidThermalData::FluidThermalData, 
    dt,tspan,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ)
    g = 9.81
    m_dot_air = fluidThermalData.m_dot_air
    m_dot_sol = fluidThermalData.m_dot_sol
    N_fin = plateFinCircularTube.N_fin
    N_row = plateFinCircularTube.N_row
    N_tube_per_row = plateFinCircularTube.N_tube_per_row
    H = plateFinCircularTube.H
    FD = plateFinCircularTube.FD
    δ_fin = plateFinCircularTube.δ_fin
    D_tube_outside = plateFinCircularTube.D_tube_outside
    Tₛₒₗ_ᵢₙ = fluidThermalData.T_sol_in
    ξₛₒₗ_ᵢₙ = fluidThermalData.X_sol_in
    IL = fluidThermalData.IL
    Tₐ_ᵢₙ = fluidThermalData.T_air
    T_wb_air = fluidThermalData.T_wb_air
    FS = plateFinCircularTube.FS
    Q = fluidThermalData.Q
    σ = plateFinCircularTube.σ
    Le = fluidThermalData.Le
    # ========================================
    ṁₐ = (m_dot_air / N_fin) * 0.5 # mass flow rate for half of the fin space
    ṁₛₒₗ_ᵢₙ = (m_dot_sol / N_fin) * 0.5 # mass flow rate for half of the fin space
    N_tube = N_tube_per_row * N_row
    H_adjuasted = (H * FD - N_tube_per_row * π * 0.25 * D_tube_outside^2) / FD
    # ========================================
    ρₛₒₗ = _ρₛₒₗ(Tₛₒₗ_ᵢₙ,ξₛₒₗ_ᵢₙ,IL)
    μₛₒₗ = _μₛₒₗ(Tₛₒₗ_ᵢₙ,ξₛₒₗ_ᵢₙ,IL)
    νₛₒₗ = μₛₒₗ / ρₛₒₗ
    δₛₒₗ = ∛(3 * ṁₛₒₗ_ᵢₙ  * νₛₒₗ / (ρₛₒₗ * g * FD))
    @show δₛₒₗ
    Uₛₒₗ_ᵣ = ṁₛₒₗ_ᵢₙ  / (ρₛₒₗ * δₛₒₗ * FD)
    ARₛₒₗ = δₛₒₗ / H_adjuasted
    Reₛₒₗ = Uₛₒₗ_ᵣ * δₛₒₗ / νₛₒₗ
    @show 4*Reₛₒₗ
    𝑘ₛₒₗ = _𝑘ₛₒₗ(Tₛₒₗ_ᵢₙ,ξₛₒₗ_ᵢₙ,IL)
    cpₛₒₗ = _cpₛₒₗ(Tₛₒₗ_ᵢₙ,ξₛₒₗ_ᵢₙ,IL)
    Prₛₒₗ = cpₛₒₗ * μₛₒₗ / 𝑘ₛₒₗ
    # ========================================
    ωₐ_ᵢₙ = HAPropsSI("W", "T", Tₐ_ᵢₙ, "P", 101325.0 , "Twb", T_wb_air)
    ρₐ = _ρₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    μₐ = _μₐ(Tₐ_ᵢₙ)
    νₐ = μₐ / ρₐ
    𝑘ₐ = _kₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    αₐ = _αₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    cpₐ = _cpₐ(Tₐ_ᵢₙ, ωₐ_ᵢₙ)
    Prₐ = νₐ / αₐ
    δₐ = 0.5FS - δₛₒₗ
    @show δₐ
    Uₐ_ᵣ = ṁₐ / (ρₐ * δₐ * FD)
    Reₐ = Uₐ_ᵣ * δₐ / νₐ
    ARₐ = δₐ / H_adjuasted
    uᵢₙₜ = 0.5g * δₛₒₗ^2 / νₛₒₗ
    dpdx = -(3.0 * μₐ * uᵢₙₜ / (δₐ^2)) - (3.0 * μₐ * ṁₐ / (ρₐ * (δₐ ^ 3) * FD))
    # ========================================
    ∂Qᵣ = (Q / N_fin) * 0.5
    @show ∂Qᵣ
    ṁₐᵢᵣ_ᵢₙ = ṁₐ
    A_c = (FS - δₛₒₗ) * FD - N_tube_per_row * (FS - δₛₒₗ) * D_tube_outside
    uₘₐₓ = ṁₐ / (ρₐ * A_c)
    @show uₘₐₓ
    iₛₒₗ_ᵢₙ = _iₛₒₗ(Tₛₒₗ_ᵢₙ,ξₛₒₗ_ᵢₙ,IL)
    iₐ_ᵢₙ = HAPropsSI("H", "T", Tₐ_ᵢₙ, "P", 101325.0 , "Twb", T_wb_air)
    NTUᴰₐᵢᵣ = NTU(ρₐ, uₘₐₓ, D_tube_outside, μₐ, Prₐ, FS,
            FD, δₛₒₗ, H_adjuasted, N_tube, N_row, 𝑘ₐ, cpₐ, Le, ṁₐ,δ_fin)
    @show NTUᴰₐᵢᵣ
    # ========================================
    solve_coil_ode!(IL ,H_adjuasted ,Le ,∂Qᵣ ,ṁₐᵢᵣ_ᵢₙ ,NTUᴰₐᵢᵣ ,σ ,ṁₛₒₗ_ᵢₙ ,ξₛₒₗ_ᵢₙ ,iₛₒₗ_ᵢₙ , ωₐ_ᵢₙ, iₐ_ᵢₙ,
                    dt,tspan,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ)
    nothing
end
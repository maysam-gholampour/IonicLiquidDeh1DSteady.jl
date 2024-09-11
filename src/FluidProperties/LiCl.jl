"""
    Manuel R. Conde,
    Properties of aqueous solutions of lithium and calcium chlorides: formulations for use in air conditioning equipment design,
    International Journal of Thermal Sciences,
    Volume 43, Issue 4,
    2004,
    Pages 367-382,
    ISSN 1290-0729,
    https://doi.org/10.1016/j.ijthermalsci.2003.09.003.
    (https://www.sciencedirect.com/science/article/pii/S1290072903001625)
    Abstract: The dehydration of air, for air conditioning purposes, either for human comfort or for industrial processes, is done most of the times by making it contact a surface at a temperature below its dew point. In this process not only is it necessary to cool that surface continuously, but also the air is cooled beyond the temperature necessary to the process, thus requiring reheating after dehumidification. Although the equipment for this purpose is standard and mostly low-cost, the running costs are high and high grade energy is dissipated at very low efficiency. Alternative sorption-based processes require only low grade energy for regeneration of the sorbent materials, thus incurring lower running costs. On the other hand, sorption technology equipment is usually more expensive than standard mechanical refrigeration equipment, which is essentially due to their too small market share. This paper reports the development of calculation models for the thermophysical properties of aqueous solutions of the chlorides of lithium and calcium, particularly suited for use as desiccants in sorption-based air conditioning equipment. This development has been undertaken in order to create consistent methods suitable for use in the industrial design of liquid desiccant-based air conditioning equipment. We have reviewed sources of measured data from 1850 onwards, and propose calculation models for the following properties of those aqueous solutions: Solubility boundary, vapour pressure, density, surface tension, dynamic viscosity, thermal conductivity, specific thermal capacity and differential enthalpy of dilution.
    Keywords: Liquid desiccants; Properties; Air conditioning; Open absorption; Lithium chloride; Calcium chloride; Calculation models
"""

function f_Pᵥₐₚₒᵣ_ₛₒₗ(ξ,θ,::LiCl)
    π₀ = 0.28
    π₁ = 4.30
    π₂ = 0.60
    π₃ = 0.21
    π₄ = 5.10
    π₅ = 0.49
    A = 2.0 - (1.0 + (ξ / π₀)^π₁)^π₂
    B = (1.0 + (ξ / π₃)^π₄)^π₅ - 1.0
    return A + B * θ
end

function _Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ,::LiCl)
    π₆ = 0.362
    π₇ = -4.75
    π₈ = -0.40
    π₉ = 0.03
    T_c_H2O = 647.226
    π₂₅ = 1.0 - (1.0 + (ξ / π₆)^π₇)^π₈ - π₉ * exp(- (ξ - 0.1)^2 / 0.005)
    θ = T / T_c_H2O
    P = f_Pᵥₐₚₒᵣ_ₛₒₗ(ξ,θ,LiCl()) * P_H2O(T) * π₂₅
    P
end

function _ρₛₒₗ(T, ξ,::LiCl)
    ρ₀ = 1.0
    ρ₁ = 0.540966
    ρ₂ = -0.303792
    ρ₃ = 0.100791
    T_c_H2O = 647.226
    θ = T / T_c_H2O
    τ = 1 - θ
    return ρ_H2O(τ) * (ρ₀ + ρ₁ * (ξ / (1 - ξ)) + ρ₂ * (ξ / (1 - ξ))^2 + ρ₃ * (ξ / (1 - ξ))^3)  
end

function _μₛₒₗ(T, ξ,::LiCl)
    ξ_ = ξ / ((1.0 - ξ) ^ (1.0 /0.6))
    η₁ = 0.090481
    η₂ = 1.390262
    η₃ = 0.675875
    η₄ = -0.583517
    T_c_H2O = 647.226
    θ = T / T_c_H2O
    return η_H2O(θ) * exp(η₁ * ξ_ ^ 3.6 + η₂ * ξ_ + η₃ * (ξ_ / θ) + η₄ * (ξ_ ^ 2))
end

function _cpₛₒₗ(T, ξ,::LiCl)
    A = 1.43980 
    B = -1.24317 
    C = -0.12070 
    D = 0.12825 
    E = 0.62934 
    F = 58.5225 
    G = -105.6343 
    H = 47.7948
    f1 = (A * ξ + B * ξ^2 + C * ξ^3) * (ξ ≤ 0.31) + (D + E * ξ) * (ξ > 0.31)
    θ = (T / 228.0) - 1.0
    f2 = F * (θ ^ 0.02) + G * (θ ^ 0.04) + H * (θ ^ 0.06) 
    Cpₛₒₗ = cp_H2O(T) * (1.0 - f1 * f2)
    Cpₛₒₗ
end

function _𝑘ₛₒₗ(T, ξ,::LiCl)
    Iₛ = 1.0
    # ================================
    # https://iopscience.iop.org/article/10.1088/1757-899X/562/1/012102
    M = 42.394
    # ================================
    α₀ = 10.8958e-3
    α₁ = -11.7882e-3
    αᵣ = α₀ + α₁ * ξ
    ξₑ = ξ * _ρₛₒₗ(T, ξ,LiCl()) * Iₛ / M
    λₛₒₗ = λ_H2O(T) - αᵣ * ξₑ
    λₛₒₗ
end 

function _iₛₒₗ(T, ξ,::LiCl)
    H₁ = 0.845 
    H₂ = -1.965 
    H₃ = -2.265 
    H₄ = 0.6 
    H₅ = 169.105 
    H₆ = 457.850 
    T_c_H2O = 647.226
    θ = T / T_c_H2O
    Δh_d0 = H₅ + H₆ * θ
    if ξ > 0.599
        ξ = 0.599
    end
    ξ = maximum(0.0,ξ)
    ξ_ = ξ / (H₄ - ξ)
    Δh_d = Δh_d0 * (1 + (ξ_ / H₁) ^ H₂) ^ H₃
    return Δh_d * 1e3
end

function _σₛₒₗ(T, ξ,::LiCl)
    T_c_H2O = 647.226
    θ = T / T_c_H2O
    σ₁ = 2.757115
    σ₂ = -12.011299  
    σ₃ = 14.751818 
    σ₄ = 2.443204 
    σ₅ = -3.147739
    σₛₒₗ = σ_H2O(θ) * (1 + 
            σ₁ * ξ +
            σ₂ * ξ * θ +
            σ₃ * ξ * θ^2 +
            σ₄ * ξ^2 +
            σ₅ * ξ^3)
    σₛₒₗ
end

function calculate_T_sol(iᵛₛₒₗ, ξ,::LiCl ;T_lower=-20.0 + 273.15, T_upper=95.0 + 273.15) 
    f(T, p)= _iₛₒₗ(T, p[2],LiCl()) - p[1]
    p = @SVector[iᵛₛₒₗ,ξ]
    T_span = @SVector[T_lower , T_upper]
    prob = IntervalNonlinearProblem(f, T_span, p)
    result = solve(prob, ITP())
    return calculate_T_barrier(result)
end



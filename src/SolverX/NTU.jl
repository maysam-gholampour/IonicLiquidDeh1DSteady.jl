export NTU

function NTU(ρ, uₘₐₓ, D_tube, μ, Pr, fin_spacing,
    fin_depth, δₛₒₗ, H, N_tube, N_row, k, cp, Le, ṁₐ,δ_fin)  
   A_fin = 2 * (H * fin_depth - N_tube * π * D_tube^2 / 4)
   A_tube = N_tube * π * D_tube * (fin_spacing - 2 * δₛₒₗ)
   Aₒ = A_fin + A_tube
   # https://doi.org/10.1115/1.2826030
   Re_D = ρ * uₘₐₓ * D_tube / μ
   Fₚ = fin_spacing + δ_fin
   j = 0.394 * Re_D^(-0.392) * ((δ_fin + δₛₒₗ) / D_tube)^(-0.0449) * N_row^(-0.0897) * (Fₚ / D_tube)^(-0.212)  
   f = 1.039 * Re_D^(-0.418) * ((δ_fin + δₛₒₗ) / D_tube)^(-0.104) * N_row^(-0.0935) * (Fₚ / D_tube)^(-0.197)  

   Nu = j * Re_D * Pr^(1/3)
   hₒ = Nu * k / D_tube
   hᴰₐᵢᵣ = hₒ / cp * Le^(2.0/3.0)
   NTUᴰₐᵢᵣ = hᴰₐᵢᵣ * fin_depth * H / ṁₐ
   NTUᴰₐᵢᵣ
   #= 
   # https://doi.org/10.1115/1.4064329
   Xₘ = Pₜ * 0.5
   Xₗ = 0.5 * √((0.5 * Pₜ)^2  + Pₗ^2)
   r = D_tube_inside * 0.5
   R_eq = 1.27 * Xₘ * √((Xₗ / Xₘ) - 0.3)
   Φ = ((R_eq / r) - 1.0) * (1.0 + 0.35 * log(R_eq / r))
   m = √(2 * hₒ / (k_fin * δ_fin))
   η = tanh(m * r * Φ) / (m * r * Φ)
   ηₒ = 1 - (A_fin / Aₒ) * (1 - η) 
   =#

end
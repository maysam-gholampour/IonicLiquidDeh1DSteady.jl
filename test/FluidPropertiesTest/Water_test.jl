@testset "Water properties" begin
    @test iᵥ_ₛₐₜ(300) ≈ 2.550049e6 atol=1.0
    @test i_fg(300) ≈ 2.437285e6 atol=1.0
end
    
@testset "Manuel R. Conde" begin
    @test P_H2O(100.0 + 273.15) ≈ 101146.0 atol=1.0

    T_c_H2O = 647.226
    T = 273.15 + 25.0
    θ = T / T_c_H2O
    τ = 1 - θ

    @test ρ_H2O(τ) ≈ 984.025 atol=1e-3
    @test η_H2O(θ) ≈ 0.00089 atol=1e-5
    @test cp_H2O(T) ≈ 4142.9 atol=0.1
    @test λ_H2O(T) ≈ 0.6065 atol=1e-4
    @test σ_H2O(θ) ≈ 0.07198 atol=1e-5
end




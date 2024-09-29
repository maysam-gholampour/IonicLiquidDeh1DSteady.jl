begin "CreCOPlus5100 Properties"
    T = 35.0 + 273.15
    ξ = 0.7
    # @test _Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 3127.884 atol=1e-3 # This test is for table interpolation
    @test _Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 3127.884 * 1.03 atol=1e-3
    @test _ρₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 1117.0205 atol=1e-4
    @test _μₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 0.0122699 atol=1e-5
    @test _cpₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 2429.7558 atol=1e-4
    @test _𝑘ₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 0.3268 atol=1e-4
    @test _Δh(T, ξ, CreCOPlus5100()) ≈ -71999.9999 atol=1e-4 # This test is for table interpolation
    # @test _Δh(T, ξ, CreCOPlus5100()) ≈ -77807.0801 atol=1e-4
    @test _σₛₒₗ(T, ξ, CreCOPlus5100()) ≈ 0.0391 atol=1e-4
    i = _iₛₒₗ(T, ξ, CreCOPlus5100())
    @test calculate_T_sol(i, ξ, CreCOPlus5100()) ≈ 35.0 + 273.15 atol=1e-3
end


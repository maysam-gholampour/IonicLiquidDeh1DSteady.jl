@testset "Water properties" begin
    @test iᵥ_ₛₐₜ(300) ≈ 2.550049e6 atol=1.0
    @test i_fg(300) ≈ 2.437285e6 atol=1.0
end
@testset "Air properties" begin
    @test _μₐ(300) ≈ 1.8537e-5 atol=1e-8
    @test _kₐ(300,0.0) ≈ 0.02638 atol=1e-4
    @test _cpₐ(300,0.0) ≈ 1006.35 atol=0.01
    @test _ρₐ(300,0.0) ≈ 1.177 atol=1e-3
    @test _νₐ(300,0.0) ≈ 1.5749e-5 atol=1e-9
    @test _αₐ(300,0.0) ≈ 2.22e-5 atol=1e-7
end